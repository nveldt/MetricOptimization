

## This works for general (dense) correlation clustering problems
function Dykstra_CD_TFA(A::SparseMatrixCSC{Float64,Int64},GapTol::Float64=1e-3,ConTol::Float64=1e-3,
                        filename::String="DykstraCDoutput",gam::Float64=10.0,maxits::Int64=1000,statusFrequency::Int64=20,stagnationTol::Float64=1e-12)
        n = size(A,1)
        open(filename, "w") do f
                write(f, "Output from Hildreth's (i.e. Dykstra's) LambdaCC Triangle Fixing Algorithm\n")
                write(f, "Gamma = $gam, Primal/Dual Gap tolerance = $GapTol, Constraint Tolerance = $ConTol \n")
        end

        # E = "error", eij = (xij - dij)
        # F = |E| i.e. fij = |eij| at optimum
        # D is the "desired distance" or just 1-Aij (anti-adjacency)
        # W is the weights matrix for LambdaCC, 1-lam or lam for each entry
        E = zeros(n,n)
        F = -gam*ones(n,n)

        # Matrices for double-loop correction terms (purposely dense)
        P = zeros(n,n)
        Q = zeros(n,n)

        # Allocate space for the output
        LPobjs = Vector{Float64}()
        duals = Vector{Float64}()
        primals = Vector{Float64}()
        gaps = Vector{Float64}()
        ConViolation = Vector{Float64}()

        # Correction term vector for triangle constraints
        current_corrections = Vector{Tuple{Int64,Float64}}()
        next_corrections = Vector{Tuple{Int64,Float64}}()
        push!(current_corrections,(0,0.0))

        # First time through constraints
        cyclic_triple_loop!(D,E,W,current_corrections,next_corrections)
        cyclic_double_loop!(E,F,P,Q,W)

        iter = 0

        # Keep track of worst constraint violations
        lastConCheck = 1.0
        lastTriCheck = 1.0

        # This keeps track of the dual variables times the right hand side.
        # It is used for computing the dual objective function.
        BtDual = 0.0

        # Initialize outside while loop
        X = zeros(n,n)
        FinalCon = 0.0
        Finalits = 0.0
        Finalobj = 0.0
        FinalGap = 0.0
        R = 0.0         # Ratio between the LP part and 2-norm penalty in objective
        Bty = 0.0

        while true
                iter += 1

                # Update which corretions to make in this iteration
                current_corrections = next_corrections
                next_corrections = Vector{Tuple{Int64,Float64}}()
                push!(current_corrections,(0,0.0))

                # Visiting the O(n^3) metric constraints is the bottleneck
                tic()
                BtDual = cyclic_triple_loop!(D,E,W,current_corrections,next_corrections)
                TriTime = toq()

                cyclic_double_loop!(E,F,P,Q,W)


                # Check convergence and give a progress update
                ConCheck, lastConCheck, lastTriCheck, objective, LPobjs, primals, duals, gaps, ConViolation,
                gap, stagnated, roundconverge, roundR, roundedGap, R, Bty =
                report_progress_CCDykstra(A,D,E,F,W,current_corrections,filename,iter,ConTol,
                        lastConCheck,lastTriCheck,BtDual,TriTime,gam,LPobjs,
                        primals,duals,gaps,ConViolation,statusFrequency,stagnationTol)

                # Stop if a desired constraint tolerance AND relative
                # primal/dual gap is small enough
                if abs(gap) < GapTol && ConCheck
                  println("Converged to within tolerance after $iter iterations")
                  X = E+D
                  Finalobj = LPobj(abs.(X-D),W)
                  Finalits = iter
                  FinalGap = gap
                  break
                end

                # If progress stagnates, return
                if stagnated
                        println("Progress stagnated, returned early")
                        X = E+D
                        Finalobj = LPobj(abs.(X-D),W)
                        FinalGap = gap
                        break
                end

                # If rounding the current iterate produces a solution, stop and
                # return that rounded solution
                if roundconverge
                        println("\t Converged to within desired tolerance by rounding the solution")
                        r = roundR
                        X = round.(E+D,r)
                        PrimalRound = Wnorm(W,X-D,abs.(X-D))/(2*gam) + LPobj(abs.(X-D),W)
                        FinalCon = FullTriangleCheck(X)
                        FinalGap = roundedGap
                        Finalobj = LPobj(abs.(X-D),W)
                        break
                end

                if iter >= maxits
                        println("Reached maximum number of iterations")
                        X = E+D
                        Finalobj = LPobj(abs.(X-D),W)
                        FinalGap = gap
                end

        end

        Finalits = iter
        return X, LPobjs, primals, duals, gaps, ConViolation, FinalCon, Finalits, Finalobj, FinalGap, R, Bty
end

# Collect variables and report progress each iteration
function report_progress_CCDykstra(A::SparseMatrixCSC{Float64,Int64},D::Matrix{Float64},E::Matrix{Float64},F::Matrix{Float64},W::Matrix{Float64},current_corrections::Vector{Tuple{Int64,Float64}},
    filename::String,iter::Int64,ConTol::Float64,lastConCheck::Float64,lastTriCheck::Float64,BtDual::Float64,triTime::Float64,
    gam::Float64,LPobjs::Vector{Float64},primals::Vector{Float64},duals::Vector{Float64},gaps::Vector{Float64},ConViolation::Vector{Float64},statusFrequency::Int64,stagnationTol::Float64)

        n = size(A,1)
        # Check whether constraint tolerance for metric constraints is satisfied
        tricheck = TriangleCheck(D+E,ConTol)

        # Constraint check for non-metric constraints
        doublecheck = DoubleCheck(E,F,ConTol)

        # True or false--are constraints satisfied to within tolerance?
        ConCheck = tricheck*doublecheck

        # Keep an updated record of the LP objective corresponding to the X variables
        #objective = LPcc_obj(A,D+E,lam)
        objective = LPobj(abs.(E),W)

        push!(LPobjs,objective)

        # Compute the objective with the F vector, which in theory
        # at the optimum will be the same, but isn't exaclty during the process
        fobj = LPobj(F,W)

        # x'*W*x where x = [e, f]
        xWx = Wnorm(W,E,F)

        # Primal objective = 1/(2*gamma)*x'*W*x + c'*x; the c'*x is the fobj function
        PrimalQP = xWx/(2*gam) + fobj
        R = xWx/(2*gam*fobj)
        # Dual objective, -b'*y - 1/(2*gamma)x'*W*x
        DualQP = -BtDual/gam - xWx/(2*gam)
        Bty = -BtDual/gam
        gap = (PrimalQP-DualQP)/DualQP

        # Round solutions for readability in output
        gapround = round(gap,5)
        PriRound = round(PrimalQP,5)
        DuRound = round(DualQP,5)
        ob = round(objective,3)
        time = round(triTime,3)

        if iter > 2
                stagnated = stagnationCheck(primals[end],PrimalQP,duals[end],DualQP,stagnationTol)
        else
                stagnated = false
        end
        # Check how dense the dual variable vector is--how many of (n choose 3)
        # triangle constraints have a nonzero correction term?
        nnz = length(current_corrections)
        ysparsity = round(nnz*6/(n*(n-1)*(n-2)),5)

        # Do a full constraint check every so often, which takes longer
        # than just answering true or false whether constraints are satisfied
        roundconverge = false
        specialPrint = false
        PrimalRound = 1.0
        roundTri = 1.0
        roundedGap = 1.0
        GapRoundTol = 1e-2
        ConRoundTol = 1e-2
        roundR = 0
        roundedGap = 0
        if iter%statusFrequency == statusFrequency-1
          lastConCheck, lastTriCheck = FullConstraintCheck(D+E,E,F)

          # If we are close to optimal, try rounding the solution vector
          # and checking feasibility and the new gap. Stop early if this
          # rounded solution is good enough.

          #tic()
          if abs(gap) < GapRoundTol && lastConCheck < ConRoundTol
                # Round to a few decimals. If any one is good enough, stop
                for r = 3:-1:1
                        Xr = round.(E+D,r)
                        xWxr = Wnorm(W,Xr-D,abs.(Xr-D))
                        objr = LPobj(abs.(Xr-D),W)
                        PrimalRound = xWxr/(2*gam) + objr
                        R = xWxr/(2*gam*objr)
                        roundTri = FullTriangleCheck(Xr)
                        #@show roundTri
                        roundedGap = (PrimalRound-DualQP)/DualQP
                        if roundTri < ConTol && abs(roundedGap) < GapTol
                                E = Xr - D
                                roundconverge = true
                                roundR = r  # save the
                                break
                        end
                end
                specialPrint = true
          end
          #timeRound = toc()
        end

        # Save current progress
        push!(primals,PrimalQP)
        push!(duals,DualQP)
        push!(gaps,gap)
        push!(ConViolation,lastConCheck)


        println("$iter:\t Dual = $DuRound, Primal = $PriRound, gap = $gapround, ConVio: $lastConCheck, TriVio = $lastTriCheck, 3Loop: $time, LPobj: $ob, Bty: $Bty")
        if specialPrint
                PR = round(PrimalRound,5); rgap = round(roundedGap,5); rTri = round(roundTri,5)
                println("\t Rounded: \tNewPrimal = $PR, gap = $rgap, ConVio = $rTri")
        end

        open(filename, "a") do f
          write(f,"$iter:\t Dual = $DuRound, Primal = $PriRound, gap = $gapround, ConVio: $lastConCheck, TriVio = $lastTriCheck, 3Loop: $time, LPobj: $ob, nnz(y) = $nnz, sparsity = $ysparsity\n")
        end

        return ConCheck, lastConCheck, lastTriCheck, objective, LPobjs, primals, duals, gaps, ConViolation, gap, stagnated, roundconverge, roundR, roundedGap, R, Bty
end


# Visits double loop constraints for LambdaCC LP
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
    delta = -Eij - Fij

    if delta > 0.0
      E[j,i] = Eij + delta/2
      F[j,i] = Fij + delta/2
      P[j,i] = delta/2
    else
      P[j,i] = 0.0
      E[j,i] = Eij
      F[j,i] = Fij
    end
end
end

for i = 1:n-1
        for j = i+1:n
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

# Visit each triangle constraint and perform necessary corrections and projections
# as a part of Dykstra's method
function cyclic_triple_loop!(D::Matrix{Float64},E::Matrix{Float64},W::Matrix{Float64},now_corrections::Vector{Tuple{Int64,Float64}},next_corrections::Vector{Tuple{Int64,Float64}})
  n = size(D,1)
  epsi = 0

  correctionsLength = length(now_corrections)
  nowInd = 1

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
