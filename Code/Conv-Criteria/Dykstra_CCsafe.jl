#
# This code is updated with new and better convergence criteria:
#       * Constraints must be satisfied to within a given AppTolerance
#       * The approximation factor is checked
#
#       Coded by Nate Veldt on March 20, 2018

include("ConvCrit_Helper.jl")

# DYKSTRA_LAMCC_TFA
#
# A Triangle Fixing Algorithm for the LambdaCC LP relaxation
#     based on Dykstra's projection algorithm
#
function Dykstra_lamCC_TFA(A::SparseMatrixCSC{Float64,Int64},AppTol::Float64=1e-3,ConTol::Float64=1e-3,
                lam::Float64=0.5,filename::String="DykstraLamCCoutput",gam::Float64=10.0,maxits::Int64=1000)

        n = size(A,1)
        open(filename, "w") do f
                write(f, "Output from Hildreth's/Dykstra's lamCC TFA\n")
                write(f, "Lambda = $lam, gamma = $gam, Approximation Tolerance = $AppTol, Constraint Tolerance = $ConTol \n")
        end

        # E = "error", eij = (xij - dij)
        # F = |E| i.e. fij = |eij| at optimum
        # D is the "desired distance" or just 1-Aij (anti-adjacency)
        # W is the weights matrix for LambdaCC, 1-lam or lam for each entry
        E = zeros(n,n)
        F = -gam*ones(n,n)
        D,W = ConstructDandW(A)

        # Matrices for double-loop correction terms (purposely dense)
        P = zeros(n,n)
        Q = zeros(n,n)

        # Correction term vector for triangle constraints
        current_corrections = Vector{Tuple{Int64,Float64}}()
        next_corrections = Vector{Tuple{Int64,Float64}}()

        push!(current_corrections,(0,0.0))

        # First time through constraints
        hildreth_triple_loop!(D,E,W,current_corrections,next_corrections)
        hildreth_double_loop!(E,F,P,Q,W)

        iter = 0

        # Variable storing the most recent check on constraint satisfaction
        lastConCheck = 1.0
        lastTriCheck = 1.0

        SlackMax = 0.0
        TotalSlack = 0.0

        # This keeps track of the dual variables times the right hand side.
        # It is used for computing the dual objective function.
        BtDual = 0.0

        # This is for the correlation clustering LP objective function,
        # which isn't technically what we're minimizing here
        objective = 0.0
        Ratio = 0.0
        while iter < maxits
                iter += 1
                vecStart = [vec(E);vec(F)]

                # An empty vector to fill future correction terms in
                current_corrections = next_corrections
                next_corrections = Vector{Tuple{Int64,Float64}}()
                push!(current_corrections,(0,0.0))

                # Visiting the O(n^3) metric constraints is the bottleneck
                tic()
                BtDual = hildreth_triple_loop!(D,E,W,current_corrections,next_corrections)
                TriTime = toq()

                # The other double loop doesn't update dual variables correctly
                hildreth_double_loop!(E,F,P,Q,W)

                # Poor test for convergence: checking change in solution iterate
                vecNext = [vec(E);vec(F)]
                vecChange = sqrt(mapreduce(z -> abs(z[1]-z[2])^2, +, zip(vecStart,vecNext)))

                # Check convergence and give a progress update
                X = E[:,:]+D[:,:]
                ConCheck, lastConCheck, lastTriCheck, objective, Ratio, SlackMax, TotalSlack = report_progress(A,D,X,E,F,P,Q,next_corrections,filename,vecChange,iter,ConTol,lam,lastConCheck,lastTriCheck,BtDual,TriTime,gam,SlackMax,TotalSlack)

                # We stop if the constraints are sufficiently satisfied
                # and the approximation ratio is small enough

                if Ratio < AppTol && ConCheck && Ratio > .5
                  break
                end

        end
        X = E+D
        FinalCon, Finaltr = FullConstraintCheck(X,E,F)
        Finalobj = LPcc_obj(A,E+D,lam)
        Finalits = iter
        return X, Finaltr, Finalobj, Finalits, Ratio
end

function report_progress(A::SparseMatrixCSC{Float64,Int64},D::Matrix{Float64},X::Matrix{Float64},E::Matrix{Float64},F::Matrix{Float64},P::Matrix{Float64},Q::Matrix{Float64},current_corrections::Vector{Tuple{Int64,Float64}},
    filename::String,xChange::Float64,iter::Int64,ConTol::Float64,lam::Float64,lastConCheck::Float64,lastTriCheck::Float64,BtDual::Float64,triTime::Float64,gam::Float64,SlackMax::Float64,TotalSlack::Float64)

        # Check to see if the constraint tolerance is satisfied for metric constraints
        # This is just true or false
        tricheck = TriangleCheck(X,ConTol)

        # Constraint check for non-metric constraints
        doublecheck = DoubleCheck(E,F,ConTol)

        # True or false--are constraints satisfied to within tolerance?
        ConCheck = tricheck*doublecheck

        # Keep an updated record of the LP objective corresponding to the X variables
        objective = LPcc_obj(A,X,lam)

        # Compute the objective with the F vector, which in theory
        # at the optimum will be the same, but isn't exaclty during the process
        fobj = LPobj(F,W)

        # x'*W*x where x = [e, f]
        xWx = Wnorm(W,E,F)

        # Primal objective = 1/(2*gamma)*x'*W*x + c'*x; the c'*x is the fobj function
        PrimalQP = round(xWx/(2*gam) + fobj,10)

        # Dual objective, -b'*y - 1/(2*gamma)x'*W*x
        DualQP = round(-BtDual/gam - xWx/(2*gam),10)

        if DualQP > 0 && PrimalQP > 0
                Ratio = round(PrimalQP/DualQP,2)
        else
                Ratio = 0
        end

        ch = round(xChange,3)
        ob = round(objective,3)
        time = round(triTime,3)

        # Do a full constraint check every so often, which takes longer
        # than just answering true or false whether constraints are satisfied
        #if iter%20 == 19
          lastConCheck, lastTriCheck = FullConstraintCheck(X,E,F)

          # If desired, we can check how well KKT conditions are satisfied
          BtDual2, SlackMax, TotalSlack, yTAxb = CheckKKT(X,D,E,F,current_corrections,P,Q)
        #end
        #@show BtDual, BtDual2, yTAxb
        nnz = length(current_corrections)
        #println("Iter $iter: ||x_k - x_{k+1}|| = $ch, 3Loop = $time, Obj: $ob, ConVio = $lastConCheck, SlackVio = $SlackMax , T-Slack = $TotalSlack")
        println("Iter $iter:\t Dual = $DualQP, Primal = $PrimalQP, ConSat: $lastConCheck, TriSat = $lastTriCheck, Ratio = $Ratio , 3Loop: $time, LPobj: $ob, nnz = $nnz, xWx = $xWx, BtDual = $BtDual, SlackVio = $SlackMax , T-Slack = $TotalSlack")
        println("$iter:\t Dual = $DualQP, Primal = $PrimalQP, ConSat: $lastConCheck, ApproxRatio = $Ratio, 3Loop: $time, LPobj: $ob")

        open(filename, "a") do f
          write(f,"$iter:\t Dual = $DualQP, Primal = $PrimalQP, ConSat: $lastConCheck, ApproxRatio = $Ratio, 3Loop: $time, LPobj: $ob\n")
        end

        return ConCheck, lastConCheck, lastTriCheck, objective, Ratio, SlackMax, TotalSlack
end


# Visits double loop constraints for LambdaCC LP. Does NOT update the
# dual variables correctly, so the primal variables will converge if you use
# this function but the Primal/Dual objective values will be off.
function cyclic_double_loop!(E::Matrix{Float64},F::Matrix{Float64},
    P::Matrix{Float64},Q::Matrix{Float64},W::Matrix{Float64})

 n = size(P,1)
 @inbounds for i = 1:n-1
  for j = i+1:n

    #  -E - F <= 0
    # corrections
    cor = P[j,i]
    Eij = E[j,i] - cor
    Fij = F[j,i] - cor
    delta = - Eij - Fij

    if delta > 0.0
      E[j,i] = Eij + delta/2
      F[j,i] = Fij + delta/2
      P[j,i] = delta/2
    else
      E[j,i] = Eij
      F[j,i] = Fij
      P[j,i] = 0.0
    end

    #  E - F <= 0
    # corrections
    cor = Q[j,i]
    Eij = E[j,i] + cor
    Fij = F[j,i] - cor
    delta = Eij - Fij
    if delta > 0.0
      E[j,i] = Eij - delta/2
      F[j,i] = Fij + delta/2
      Q[j,i] = delta/2
    else
      Q[j,i] = 0.0
      E[j,i] = Eij
      F[j,i] = Fij
    end

  end
end

end

# Hildreth style double loop. In theory does the exact same thing as
# the cyclic_double_loop function, but this one updates the dual variables
# correctly.
function hildreth_double_loop!(E::Matrix{Float64},F::Matrix{Float64},
    P::Matrix{Float64},Q::Matrix{Float64},W::Matrix{Float64})

 n = size(P,1)
 @inbounds for i = 1:n-1
  for j = i+1:n

    #  -E - F <= 0
    thetaI = (-E[j,i] - F[j,i])*W[j,i]/2

    delta = min(-thetaI,P[j,i])
    E[j,i] = E[j,i] - delta/W[j,i]
    F[j,i] = F[j,i] - delta/W[j,i]
    P[j,i] = P[j,i] - delta

    #  E - F <= 0
    # corrections
   thetaI = (E[j,i] - F[j,i])*W[j,i]/2
   delta = min(-thetaI,Q[j,i])

   E[j,i] = E[j,i] + delta/W[j,i]
   F[j,i] = F[j,i] - delta/W[j,i]
   Q[j,i] = Q[j,i] - delta
  end
end

end

function cyclic_triple_loop!(D::Matrix{Float64},E::Matrix{Float64},W::Matrix{Float64},now_corrections::Vector{Tuple{Int64,Float64}},next_corrections::Vector{Tuple{Int64,Float64}})
  n = size(D,1)
  epsi = 0

  correctionsLength = length(now_corrections)
  nowInd = 1;
  # Grab next triplet in the list
  nextKey = now_corrections[nowInd]

  # Keep updated the dot product between dual variables and the rhs of Ax <= b
  BtDual = 0.0
  @inbounds for i = 1:n-2
    for  j = i+1:n-1
     Dij = D[j,i]
     Wij = W[j,i]
    for  k = j+1:n

    Dik = D[k,i]
    Djk = D[k,j]
    Eij = E[j,i]
    Eik = E[k,i]
    Ejk = E[k,j]
    Wik = W[k,i]
    Wjk = W[k,j]

    ### Check triangle i,j,k
    ijkKey = (i-1)*n^2+(j-1)*n+k
    # First see if this is the next triangle with a nonzero correction variable
    if ijkKey == nextKey[1]

            cor = nextKey[2]
      # We need to scale the correction since we're projecting into the minimum
      # weighted 2-norm direction

      denom = Wij*Wjk + Wik*Wij + Wjk*Wik

      E[j,i] = Eij + cor*(Wik*Wjk/denom)
      E[k,i] = Eik - cor*(Wij*Wjk/denom)
      E[k,j] = Ejk - cor*(Wik*Wij/denom)
      Eij = E[j,i]
      Eik = E[k,i]
      Ejk = E[k,j]
      # Move along in the list of triplets with corrections
      if nowInd < correctionsLength
              nowInd +=1
              nextKey = now_corrections[nowInd]
      end
    end

    b = Dik + Djk - Dij
    mu = (Eij - Ejk - Eik - b)

    if mu > epsi
      denom = Wij*Wjk + Wik*Wij + Wjk*Wik
      E[j,i] = Eij - mu*(Wik*Wjk/denom)
      E[k,i] = Eik + mu*(Wij*Wjk/denom)
      E[k,j] = Ejk + mu*(Wik*Wij/denom)
      # Next time we see this triple we have to correct
      push!(next_corrections,(ijkKey,mu))
      BtDual += mu*b*(Wij*Wik*Wjk/denom)
    end

    ### Done checking triangle i,j,k


    ### Check triangle i,k,j
    ijkKey = (i-1)*n^2+(k-1)*n+j
    # First see if this is the next triangle with a nonzero correction variable
    if ijkKey == nextKey[1]

            cor = nextKey[2]

      denom = Wij*Wjk + Wik*Wij + Wjk*Wik
      E[j,i] = E[j,i] - cor*(Wik*Wjk/denom)
      E[k,i] = E[k,i] + cor*(Wij*Wjk/denom)
      E[k,j] = E[k,j] - cor*(Wik*Wij/denom)
      Eij = E[j,i]
      Eik = E[k,i]
      Ejk = E[k,j]
      # Move along in the list of triplets with corrections
      if nowInd < correctionsLength
              nowInd +=1
              nextKey = now_corrections[nowInd]
      end
    else
      Eij = E[j,i]
      Eik = E[k,i]
      Ejk = E[k,j]
    end

    b = -Dik + Djk + Dij
    mu = (-Eij - Ejk + Eik - b)
    if mu > epsi

      denom = Wij*Wjk + Wik*Wij + Wjk*Wik

      E[j,i] = Eij + mu*(Wik*Wjk/denom)
      E[k,i] = Eik - mu*(Wij*Wjk/denom)
      E[k,j] = Ejk + mu*(Wik*Wij/denom)

      # Next time we see this triple we have to correct
      push!(next_corrections,(ijkKey,mu))
      BtDual += mu*b*(Wij*Wik*Wjk/denom)
    end
    ### Done checking triangle i,k,j

    ### Triangle j,k,i
    ijkKey = (j-1)*n^2+(k-1)*n+i
    # First see if this is the next triangle with a nonzero correction variable
    if ijkKey == nextKey[1]

            cor = nextKey[2]

      denom = Wij*Wjk + Wik*Wij + Wjk*Wik

      E[j,i] = E[j,i] - cor*(Wik*Wjk/denom)
      E[k,i] = E[k,i] - cor*(Wij*Wjk/denom)
      E[k,j] = E[k,j] + cor*(Wik*Wij/denom)
      Eij = E[j,i]
      Eik = E[k,i]
      Ejk = E[k,j]
      # Move along in the list of triplets with corrections
      if nowInd < correctionsLength
              nowInd +=1
              nextKey = now_corrections[nowInd]
      end
    else
      Eij = E[j,i]
      Eik = E[k,i]
      Ejk = E[k,j]
    end

    b = Dik - Djk + Dij
    mu = (-Eij + Ejk - Eik - b)

    if mu > epsi
      denom = Wij*Wjk + Wik*Wij + Wjk*Wik

      E[j,i] = Eij + mu*(Wik*Wjk/denom)
      E[k,i] = Eik + mu*(Wij*Wjk/denom)
      E[k,j] = Ejk - mu*(Wik*Wij/denom)

      # Next time we see this triple we have to correct
      push!(next_corrections,(ijkKey,mu))
      BtDual += mu*b*(Wij*Wik*Wjk/denom)
    end
    ### Done checking triangle j,k,i

  end
  end
end # end triple for loop
        return BtDual
end # end function

# Hildreth-style update procedure for the double loop. Does exactly the same
# thing as the Dykstra-styple update as long as you are careful about
# updating dual variables correctly (which we do here).
function hildreth_triple_loop2!(D::Matrix{Float64},E::Matrix{Float64},W::Matrix{Float64},now_corrections::Vector{Tuple{Int64,Float64}},next_corrections::Vector{Tuple{Int64,Float64}})
  n = size(D,1)
  correctionsLength = length(now_corrections)
  nowInd = 1
  # Grab next triplet in the list
  nextKey = now_corrections[nowInd]

  BtDual = 0.0
  @inbounds for i = 1:n-2
    for  j = i+1:n-1
     Eij = E[j,i]
     Dij = D[j,i]
     Wij = W[j,i]
    for  k = j+1:n

    Dik = D[k,i]
    Djk = D[k,j]
    Eik = E[k,i]
    Ejk = E[k,j]
    Wik = W[k,i]
    Wjk = W[k,j]
    denom = Wij*Wjk + Wik*Wij + Wjk*Wik
    aWinva = denom/(Wjk*Wik*Wij)
    ### Check triangle i,j,k

    ijkKey = (i-1)*n^2+(j-1)*n+k
    # First see if this is the next triangle with a nonzero correction variable

    if ijkKey == nextKey[1]

            yI = nextKey[2]
      # Move along in the list of triplets with corrections
      if nowInd < correctionsLength
              nowInd +=1
              nextKey = now_corrections[nowInd]
      end
    else
        yI = 0.0
    end

    b = Dik + Djk - Dij
    thetaI = (Eij - Ejk - Eik - b)/aWinva

    delta = min(-thetaI,yI)

    if delta != 0
      E[j,i] = Eij + delta/Wij
      E[k,i] = Eik - delta/Wik
      E[k,j] = Ejk - delta/Wjk
    end

    if yI - delta > 0
        push!(next_corrections,(ijkKey,yI-delta))
        BtDual += b*(yI-delta)
    end

    ### Done checking triangle i,j,k

    Eij = E[j,i]
    Eik = E[k,i]
    Ejk = E[k,j]

    ### Check triangle i,k,j
    ijkKey = (i-1)*n^2+(k-1)*n+j
    if ijkKey == nextKey[1]

            yI = nextKey[2]
      # Move along in the list of triplets with corrections
      if nowInd < correctionsLength
              nowInd +=1
              nextKey = now_corrections[nowInd]
      end
    else
        yI = 0.0
    end

    b = -Dik + Djk + Dij
    thetaI = (-Eij - Ejk + Eik - b)/aWinva

    delta = min(-thetaI,yI)
    #delta = max(thetaI,-yI)

    if delta != 0
      E[j,i] = Eij - delta/Wij
      E[k,i] = Eik + delta/Wik
      E[k,j] = Ejk - delta/Wjk
    end

    # Update dual variable if it's not zero
    if yI - delta > 0
        push!(next_corrections,(ijkKey,yI-delta))
        BtDual += b*(yI-delta)
    end
    ### Done checking triangle i,k,j

    Eij = E[j,i]
    Eik = E[k,i]
    Ejk = E[k,j]
    ### Triangle j,k,i
    ijkKey = (j-1)*n^2+(k-1)*n+i
    if ijkKey == nextKey[1]

            yI = nextKey[2]
      # Move along in the list of triplets with corrections
      if nowInd < correctionsLength
              nowInd +=1
              nextKey = now_corrections[nowInd]
      end
    else
        yI = 0.0
    end

    b = Dik - Djk + Dij
    thetaI = (-Eij + Ejk - Eik - b)/aWinva

    delta = min(-thetaI,yI)
    #delta = max(thetaI,-yI)

    if delta != 0
      E[j,i] = Eij - delta/Wij
      E[k,i] = Eik - delta/Wik
      E[k,j] = Ejk + delta/Wjk
      # Next time we see this triple we have to correct
    end

    # Update dual variable if it's not zero
    if yI - delta > 0
        push!(next_corrections,(ijkKey,yI-delta))
        BtDual += b*(yI - delta)
    end
    ### Done checking triangle j,k,i

end
end
end # end triple for loop

        return BtDual
end # end TripleLoop! function


# This is the way that [Dax 2003} presents the Hildreth-like method,
# which he doesn't appropriately attribute to Hildreth.
function hildreth_triple_loop!(D::Matrix{Float64},E::Matrix{Float64},W::Matrix{Float64},
        now_corrections::Vector{Tuple{Int64,Float64}},next_corrections::Vector{Tuple{Int64,Float64}})
  n = size(D,1)
  correctionsLength = length(now_corrections)
  nowInd = 1
  # Grab next triplet in the list
  nextKey = now_corrections[nowInd]

  BtDual = 0.0
  @inbounds for i = 1:n-2
    for  j = i+1:n-1

     Dij = D[j,i]
     Wij = W[j,i]
    for  k = j+1:n

    Eij = E[j,i]
    Dik = D[k,i]
    Djk = D[k,j]
    Eik = E[k,i]
    Ejk = E[k,j]
    Wik = W[k,i]
    Wjk = W[k,j]
    denom = Wij*Wjk + Wik*Wij + Wjk*Wik
    aWinva = denom/(Wjk*Wik*Wij)
    ### Check triangle i,j,k

    ijkKey = (i-1)*n^2+(j-1)*n+k
    # First see if this is the next triangle with a nonzero correction variable

    if ijkKey == nextKey[1]

            yI = nextKey[2]
      # Move along in the list of triplets with corrections
      if nowInd < correctionsLength
              nowInd +=1
              nextKey = now_corrections[nowInd]
      end
    else
        yI = 0.0
    end

    b = Dik + Djk - Dij
    thetaI = (Eij - Ejk - Eik - b)/aWinva

    #delta = min(-thetaI,yI)
    delta = max(thetaI,-yI)

    if delta != 0
      E[j,i] = Eij - delta/Wij
      E[k,i] = Eik + delta/Wik
      E[k,j] = Ejk + delta/Wjk
    end

    if yI + delta != 0
        push!(next_corrections,(ijkKey,yI+delta))
        BtDual += b*(yI+delta)
    end

    ### Done checking triangle i,j,k

    Eij = E[j,i]
    Eik = E[k,i]
    Ejk = E[k,j]

    ### Check triangle i,k,j
    ijkKey = (i-1)*n^2+(k-1)*n+j
    if ijkKey == nextKey[1]

            yI = nextKey[2]
      # Move along in the list of triplets with corrections
      if nowInd < correctionsLength
              nowInd +=1
              nextKey = now_corrections[nowInd]
      end
    else
        yI = 0.0
    end

    b = -Dik + Djk + Dij
    thetaI = (-Eij - Ejk + Eik - b)/aWinva

    #delta = min(-thetaI,yI)
    delta = max(thetaI,-yI)

    if delta != 0
      E[j,i] = Eij + delta/Wij
      E[k,i] = Eik - delta/Wik
      E[k,j] = Ejk + delta/Wjk
    end

    # Update dual variable if it's not zero
    if yI + delta != 0
        push!(next_corrections,(ijkKey,yI+delta))
        BtDual += b*(yI+delta)
    end
    ### Done checking triangle i,k,j

    Eij = E[j,i]
    Eik = E[k,i]
    Ejk = E[k,j]
    ### Triangle j,k,i
    ijkKey = (j-1)*n^2+(k-1)*n+i
    if ijkKey == nextKey[1]

            yI = nextKey[2]
      # Move along in the list of triplets with corrections
      if nowInd < correctionsLength
              nowInd +=1
              nextKey = now_corrections[nowInd]
      end
    else
        yI = 0.0
    end

    b = Dik - Djk + Dij
    thetaI = (-Eij + Ejk - Eik - b)/aWinva

    #delta = min(-thetaI,yI)
    delta = max(thetaI,-yI)

    if delta != 0
      E[j,i] = Eij + delta/Wij
      E[k,i] = Eik + delta/Wik
      E[k,j] = Ejk - delta/Wjk
      # Next time we see this triple we have to correct
    end

    # Update dual variable if it's not zero
    if yI + delta != 0
        push!(next_corrections,(ijkKey,yI+delta))
        BtDual += b*(yI+delta)
    end
    ### Done checking triangle j,k,i

end
end
end # end triple for loop

        return BtDual
end # end TripleLoop! function


# Check the worst violation in the KKT conditions
function CheckKKT(X::Matrix{Float64},D::Matrix{Float64},E::Matrix{Float64},F::Matrix{Float64},now_corrections::Vector{Tuple{Int64,Float64}},P::Matrix{Float64},Q::Matrix{Float64})

    n = size(X,1)
    correctionsLength = length(now_corrections)
    nowInd = 1
    # Grab next triplet in the list
    nextKey = now_corrections[nowInd]

    maxi = 0.0
    SlackMax = 0.0
    TotalSlack = 0.0
    BtDual = 0.0
    yTAxb = 0.0
    @inbounds for i = 1:n-2
      for j = i+1:n-1
        a = X[j,i]
        dij = D[j,i]
        for k = j+1:n
          b = X[k,i]
          c = X[k,j]
          dik = D[k,i]
          djk = D[k,j]

          dijk = a - b - c
          dikj = b - c - a
          djki = c - a - b

          bijk = -dij+dik+djk
          bikj = dij-dik+djk
          bjki = dij+dik-djk

          # Check for complementary slackness

          ## Check triangle i,j,k
         ijkKey = (i-1)*n^2+(j-1)*n+k
         if ijkKey == nextKey[1]

            cor = nextKey[2]
            slackness = abs(cor*dijk)
            yTAxb +=cor*dijk
            BtDual += cor*bijk
            TotalSlack += slackness
            if slackness > SlackMax
                SlackMax = slackness
            end
            if nowInd < correctionsLength
                    nowInd +=1
                    nextKey = now_corrections[nowInd]
            end
         end

         ## Check triangle i,j,k
        ijkKey = (i-1)*n^2+(k-1)*n+j
        if ijkKey == nextKey[1]

           cor = nextKey[2]
           slackness = abs(cor*dikj)
           yTAxb +=cor*dikj
           BtDual += cor*bikj
           TotalSlack += slackness
           if slackness > SlackMax
               SlackMax = slackness
           end
           if nowInd < correctionsLength
                   nowInd +=1
                   nextKey = now_corrections[nowInd]
           end
        end

        ## Check triangle i,j,k
       ijkKey = (j-1)*n^2+(k-1)*n+i
       if ijkKey == nextKey[1]

          cor = nextKey[2]
          slackness = abs(cor*djki)
          yTAxb +=cor*djki
          BtDual += cor*bjki
          TotalSlack += slackness
          if slackness > SlackMax
              SlackMax = slackness
          end
          if nowInd < correctionsLength
                  nowInd +=1
                  nextKey = now_corrections[nowInd]
          end
       end

          # Check for worst constraint violation
        if dijk > maxi
            maxi = dijk
        elseif dikj > maxi
            maxi = dikj
        elseif djki > maxi
            maxi = djki
        end

        end
      end
    end
    #
    for i = 1:n-1
        for j = i+1:n
            Eij = E[j,i]
            Fij = F[j,i]
            if Eij - Fij > maxi
                maxi = Eij - Fij
            elseif -Eij - Fij > maxi
                maxi = -Eij - Fij
            end
            slackness = abs((Eij - Fij)*Q[j,i])
            yTAxb += (Eij - Fij)*Q[j,i]
            TotalSlack += slackness
            if slackness > SlackMax
                SlackMax = slackness
            end
            slackness = abs((-Eij - Fij)*P[j,i])
            yTAxb += (-Eij - Fij)*P[j,i]
            TotalSlack += slackness
            if slackness > SlackMax
                SlackMax = slackness
            end
        end
    end
    return BtDual, round(SlackMax,10), round(TotalSlack,10), yTAxb
  end
