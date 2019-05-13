# This file includes Dykstra-based Projection algorithms
# for Correlation clustering and Sparsest Cut.
#
# This version of the code is specifically designed for Julia 1.0
#
# Last updated Nate Veldt on January 24, 2019

using SparseArrays
using LinearAlgebra

# DYKSTRA_LAMCC_TFA
#
# A Triangle Fixing Algorithm for the LambdaCC LP relaxation
#     based on Dykstra's projection algorithm
#
# Inputs:
#
#       A = the unweighted, symmetric adjacency matrix for an undirected graph
#       GapTol = the tolerance for the relative gap between
#                primal and dual functions. Getting converngece to this tolerance
#               means returning a solution that is within factor (1+GapTol)
#               of the optimal solution
#       ConTol = constraint tolerance
#       lam = the LambdaCC resolution parameter
#       filename = where we save the output from running this algorithm
#       gam = the parameter controlling the relationship between the LP and QP
#       maxits = maximum number of iterations
#       roundFlag = true or false depeneding on whether you want to run the
#                       entrywise rounding procedure
#       rlow, rhigh = lower and upper integer--range of values to try the rounding procedure on
#       pretol = the initial tolerance that must be satisfied in order to apply the rounding procedure
#
# Outputs:
#
#       X = the distance matrix between nodes
#       LPobjs = a vector containing the LPobjective scores at each iteration
#       duals = a vector containing QP dual objective scores (lower bounds)
#       primals = a vector containing QP primal objective scores (not all upper bounds)
#       ConViolation = a vector of the max constraint violation at each iteration
function Dykstra_lamCC_TFA(A::SparseMatrixCSC{Float64,Int64},lam::Float64=0.5,
                GapTol::Float64=1e-3,ConTol::Float64=1e-3,gam::Float64=10.0,
                maxits::Int64=5000,statusFrequency::Int64=20,
                filename::String="DykstraLamCCoutput",
                roundFlag::Bool=false, rlow::Int64=1, rhigh::Int64=1, pretol::Float64=1e-10)

        D,W = LamCC_DandW(A,lam)
        Dykstra_CC_TFA(A,W,D,GapTol,ConTol,filename,gam,maxits,statusFrequency,roundFlag,rlow,rhigh,pretol)

end

function Dykstra_lamCC_dw(A::SparseMatrixCSC{Float64,Int64},lam::Float64=2.0,
        GapTol::Float64=1e-2,ConTol::Float64=1e-2,gam::Float64=10.0,
        maxits::Int64=5000,statusFrequency::Int64=20,filename::String="DykstraLamCCoutput",
        roundFlag::Bool=false, rlow::Int64=1, rhigh::Int64=1, pretol::Float64=1e-10)


        if lam > 1
            volA = sum(nonzeros(A))
            lam = 1/(volA)
        end
        D,W = LamCC_degweighted(A,lam)
        @assert(minimum(W) >= 0)
        Dykstra_CC_TFA(A,W,D,GapTol,ConTol,filename,gam,maxits,statusFrequency)

end

## This works for general (dense) correlation clustering problems
function Dykstra_CC_TFA(A::SparseMatrixCSC{Float64,Int64},W::Matrix{Float64},D::Matrix{Float64},GapTol::Float64=1e-5,ConTol::Float64=1e-5,
                        filename::String="DykstraCCoutput",gam::Float64=10.0,maxits::Int64=1000,statusFrequency::Int64=20,
                        roundFlag::Bool=false, rlow::Int64=1, rhigh::Int64=1, pretol::Float64=1e-10)
        n = size(A,1)
        open(filename, "w") do f
                write(f, "Output from DykstraCC \n")
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
        cyclic_triple_loopE!(D,E,W,current_corrections,next_corrections)
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
        R = 0.0        # Ratio between the LP part and 2-norm penalty in objective
        Bty = 0.0

        while true
                iter += 1

                # Update which corretions to make in this iteration
                current_corrections = next_corrections
                next_corrections = Vector{Tuple{Int64,Float64}}()
                push!(current_corrections,(0,0.0))

                # Visiting the O(n^3) metric constraints is the bottleneck
                start = time()
                BtDual = cyclic_triple_loopE!(D,E,W,current_corrections,next_corrections)
                TriTime = time()-start

                cyclic_double_loop!(E,F,P,Q,W)

                # Check convergence and give a progress update
                ConCheck, lastConCheck, lastTriCheck, objective, LPobjs, primals, duals, gaps, ConViolation,
                gap, roundconverge, roundR, roundedGap, R, Bty =
                report_progress_CCDykstra(A,D,E,F,W,current_corrections,filename,iter,ConTol,GapTol,
                        lastConCheck,lastTriCheck,BtDual,TriTime,gam,LPobjs,
                        primals,duals,gaps,ConViolation,statusFrequency,
                        roundFlag,rlow,rhigh,pretol)

                # Stop if a desired constraint tolerance AND relative
                # primal/dual gap is small enough
                if abs(gap) < GapTol && ConCheck
                  println("Converged to within tolerance after $iter iterations")
                  X = E+D
                  Finalobj = LPobj(abs.(X-D),W)
                  Finalits = iter
                  FinalGap = gap
                  FinalCon = FullTriangleCheck(X)
                  break
                end

                # If rounding the current iterate produces a solution, stop and
                # return that rounded solution
                if roundconverge
                  println("\t Converged to within desired tolerance by rounding the solution")
                  open(filename, "w") do f
                          write(f, "\t Converged to within desired tolerance by rounding the solution\n")
                  end
                  r = roundR
                  X = round.(E+D,digits = r)
                  PrimalRound = Wnorm(W,X-D,abs.(X-D))/(2*gam) + LPobj(abs.(X-D),W)
                  FinalCon = FullTriangleCheck(X)
                  FinalGap = roundedGap
                  Finalobj = LPobj(abs.(X-D),W)
                  break
                end

                if iter >= maxits
                  println("Reached maximum number of iterations")
                  open(filename, "w") do f
                          write(f, "\t Reached maximum number of iterations\n")
                  end
                  X = E+D
                  Finalobj = LPobj(abs.(X-D),W)
                  FinalGap = gap
                  break
                end

        end

        Finalits = iter
        return X, LPobjs, primals, duals, gaps, ConViolation, FinalCon, Finalits, Finalobj, FinalGap, R, Bty, next_corrections
end

# Collect variables and report progress each iteration
function report_progress_CCDykstra(A::SparseMatrixCSC{Float64,Int64},
  D::Matrix{Float64},E::Matrix{Float64},F::Matrix{Float64},W::Matrix{Float64},
  current_corrections::Vector{Tuple{Int64,Float64}},filename::String,
  iter::Int64,ConTol::Float64,GapTol::Float64,lastConCheck::Float64,
  lastTriCheck::Float64,BtDual::Float64,triTime::Float64,gam::Float64,
  LPobjs::Vector{Float64},primals::Vector{Float64},duals::Vector{Float64},
  gaps::Vector{Float64},ConViolation::Vector{Float64},statusFrequency::Int64,
  roundFlag::Bool, rlow::Int64, rhigh::Int64, pretol::Float64)

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
        gapround = round(gap,digits = 5)
        PriRound = round(PrimalQP,digits = 5)
        DuRound = round(DualQP,digits = 5)
        ob = round(objective,digits = 3)
        time = round(triTime,digits = 3)

        # Check how dense the dual variable vector is--how many of (n choose 3)
        # triangle constraints have a nonzero correction term?
        nnz = length(current_corrections)
        ysparsity = round(nnz*6/(n*(n-1)*(n-2)),digits = 5)

        # Do a full constraint check every so often, which takes longer
        # than just answering true or false whether constraints are satisfied
        roundconverge = false
        specialPrint = false
        PrimalRound = 1.0
        roundTri = 1.0
        roundedGap = 1.0
        GapRoundTol = 5e-1
        ConRoundTol = 5e-1
        roundR = 0
        roundedGap = 0

        if iter%statusFrequency == statusFrequency-1
          lastConCheck, lastTriCheck = FullConstraintCheck(D+E,E,F)

          # If we are close to optimal, try rounding the solution vector
          # and checking feasibility and the new gap. Stop early if this
          # rounded solution is good enough.

          if roundFlag && abs(gap) < pretol && lastConCheck < pretol
            # Round to a few decimals. If any one is good enough, stop
            for r = rhigh:-1:rlow
              Xr = round.(E+D,digits = r)
              xWxr = Wnorm(W,Xr-D,abs.(Xr-D))
              objr = LPobj(abs.(Xr-D),W)
              PrimalRound = xWxr/(2*gam) + objr
              R = xWxr/(2*gam*objr)
              roundTri = FullTriangleCheck(Xr)
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
        end

        # Save current progress
        push!(primals,PrimalQP)
        push!(duals,DualQP)
        push!(gaps,gap)
        push!(ConViolation,lastConCheck)


        println("$iter:\t Dual = $DuRound, Primal = $PriRound, gap = $gapround, ConVio: $lastConCheck, TriVio = $lastTriCheck, 3Loop: $time, LPobj: $ob")
        if specialPrint
                PR = round(PrimalRound,digits = 5); rgap = round(roundedGap,digits = 5); rTri = round(roundTri,digits = 5)
                println("\t Rounded: \tNewPrimal = $PR, gap = $rgap, ConVio = $rTri")
        end

        open(filename, "a") do f
          write(f,"$iter:\t Dual = $DuRound, Primal = $PriRound, gap = $gapround, ConVio: $lastConCheck, TriVio = $lastTriCheck, 3Loop: $time, LPobj: $ob, nnz(y) = $nnz, sparsity = $ysparsity\n")
        end

        return ConCheck, lastConCheck, lastTriCheck, objective, LPobjs, primals, duals, gaps, ConViolation, gap, roundconverge, roundR, roundedGap, R, Bty
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
function cyclic_triple_loopE!(D::Matrix{Float64},E::Matrix{Float64},W::Matrix{Float64},
  now_corrections::Vector{Tuple{Int64,Float64}},
  next_corrections::Vector{Tuple{Int64,Float64}})

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

end # end cyclict_triple_loop function

# DYKSTRA_SC
#
# Dykstra-based projection method for solving a quadratic program which is the
# regularized version of the Leighton-Rao Linear Programming relaxation for
# sparsest cut.
#
# Paramters:
#
# A = adjacency matrix of an undirected, unweighted graph
# GapTol = the desired relative duality gap tolerance we wish to achieve for convergence
# ConTol = tolerance for constraint violations
# lam, gam = paramters controlling the relationship between the original LP
#               and the quadratic program which is solved in practice here.
# statusFreqcuency = controls how often we perform a full convergence check,
#                       which involves a full check for the maximum constraint violations
#                       and includes the "entrywise rounding step" (see paper)
# maxits = maximum number of iterations to run before terminating
function DykstraSC(A::SparseMatrixCSC{Float64,Int64},GapTol::Float64=1e-3,ConTol::Float64=1e-5,
                lam::Float64=0.1,filename::String="DykstraLeightonRaoOutput",gam::Float64=10.0,
                maxits::Int64=1000,statusFrequency::Int64=10,
                roundFlag::Bool=true, rlow::Int64=2, rhigh::Int64=6, pretol::Float64=5e-1)

        n = size(A,1)
        open(filename, "w") do f
                write(f, "DykstraSC\n")
                write(f, "Lambda = $lam, gamma = $gam, tol = $GapTol, ConTol = $ConTol \n")
        end

        # Initialize X = -gam*A,
        # Wij = 1 if A_ij = 1
        # Wij = lam is A_ij = 0
        X,W = LeightonRaoQP_Initialize(A,lam,gam)

        # Matrices for double-loop correction terms (purposely dense)
        P = zeros(n,n)

        # Correction variable for constraint sum_ij xij = n
        SumCorrection = 0.0

        # Constant used in performing projections
        InvSumWeights = 0.0
        for i = 1:n-1
            for j = i+1:n
                InvSumWeights += 1/W[j,i]
            end
        end

        # Allocate space for the output
        LPobjs = Vector{Float64}()
        duals = Vector{Float64}()
        primals = Vector{Float64}()
        gaps = Vector{Float64}()
        ConViolation = Vector{Float64}()

        # Correction term vector (i.e. dual varialbes) for triangle constraints
        current_corrections = Vector{Tuple{Int64,Float64}}()
        next_corrections = Vector{Tuple{Int64,Float64}}()
        push!(current_corrections,(0,0.0))

        # First time through constraints
        cyclic_triple_loopX!(X,W,current_corrections,next_corrections)
        SumCorrection = the_sum_constraint!(X,W,SumCorrection,InvSumWeights)
        box_constraints!(X,W,P)

        iter = 0
        lastConCheck = 1.0

        # Make these variables global so they exist outside the while loop
        objective = 0.0
        R = 0.0
        FinalGap = 0.0
        FinalCon = 0.0
        Finalobj = 0.0
        Finalits = 0.0
        Bty = 0.0

        while true
          iter += 1

          # An empty vector to fill future correction terms in
          current_corrections = next_corrections
          next_corrections = Vector{Tuple{Int64,Float64}}()
          push!(current_corrections,(0,0.0))

          start = time()
          cyclic_triple_loopX!(X,W,current_corrections,next_corrections)
          SumCorrection = the_sum_constraint!(X,W,SumCorrection,InvSumWeights)
          box_constraints!(X,W,P)
          TriTime = time()-start


          tricheck, lastConCheck, objective, LPobjs, primals, duals, gaps,
          ConViolation, gap, roundconverge, roundR, R, roundedGap, roundTri, Bty = report_progress_DykstraLR(A,X,W,P,current_corrections,filename,iter,
          ConTol,GapTol,lam,lastConCheck,TriTime,SumCorrection,gam,LPobjs,primals,duals,gaps,ConViolation,statusFrequency,
          roundFlag, rlow, rhigh, pretol)

          # Return if the gap is less than tolerance and constraints
          # are satisfied to within a given tolerance
          if abs(gap) < GapTol && tricheck
            println("Converged without rounding procedure.")
            open(filename,"a") do f
               write(f, "Converged without rounding procedure.\n")
            end
            FinalCon = FullTriangleCheck(X)
            Finalobj = LR_obj(A,X)
            Finalits = iter
            FinalGap = gap
            break
          end

          # If rounding the current iterate produces a solution, stop and
          # return that rounded solution
          if roundconverge
            println("\t Converged to within desired tolerance by rounding the solution to $roundR decimals")
            open(filename,"a") do f
                write(f, "\t Converged to within desired tolerance by rounding the solution to $roundR decimals\n")
            end
            Xr = round.(X,digits = roundR)
            adjust = n/sum(tril(Xr))
            Xr = adjust*Xr
            X = Xr
            FinalCon = FullTriangleCheck(X)
            Finalobj = LR_obj(A,X)
            Finalits = iter
            FinalGap = roundedGap
            break
          end

          if iter >= maxits
            println("Maximum number of iterations reached")
            open(filename,"a") do f
                write(f, "Maximum number of iterations reached\n")
            end
            FinalCon = FullTriangleCheck(X)
            Finalobj = LR_obj(A,X)
            Finalits = iter
            FinalGap = gap
            break
          end
        end

        # Output final statistics and iterates
        return X, FinalCon, FinalGap, Finalobj, Finalits, R, LPobjs, duals, primals, gaps, ConViolation, Bty
end


function report_progress_DykstraLR(A::SparseMatrixCSC{Float64,Int64},X::Matrix{Float64},W::Matrix{Float64},P::Matrix{Float64},current_corrections::Vector{Tuple{Int64,Float64}},
    filename::String,iter::Int64,ConTol::Float64,GapTol::Float64,lam::Float64,lastConCheck::Float64,triTime::Float64,SumCorrection::Float64,gam::Float64,
    LPobjs::Vector{Float64},primals::Vector{Float64},duals::Vector{Float64},gaps::Vector{Float64},ConViolation::Vector{Float64},statusFrequency::Int64,
    roundFlag::Bool, rlow::Int64, rhigh::Int64, pretol::Float64)

        n = size(X,1)

        # True or false convergence check
        tricheck = TriangleCheck(X,ConTol)

        sumcheck = abs(n-sum(tril(X)))

        tricheck = tricheck*(sumcheck<ConTol)

        # nonzeros in dual vector
        nnzDelta = length(current_corrections)

        # Compute primal and dual objectives
        objective = LR_obj(A,X)
        xWx = xWnorm(W,X)
        BtDual = (n*SumCorrection)
        PrimalQP = xWx/(2*gam) + objective
        R = xWx/(2*gam*objective)
        DualQP = -BtDual/(gam) - xWx/(2*gam)
        gap = (PrimalQP-DualQP)/DualQP

        # Every "statusFrequency" number of iterations, fully check constraints
        # to find the magnitude of the worst violation.
        # Also, try rounding the current solution if we are close to convergence
        roundconverge = false
        specialPrint = false
        PrimalRound = 1.0
        roundTri = 1.0
        roundedGap = 1.0
        roundR = 0

        if iter%statusFrequency == statusFrequency-1
          lastConCheck = max(sumcheck,FullTriangleCheck(X))

          if roundFlag && abs(gap) < pretol && lastConCheck < pretol
            # Round to a few decimals. If any one is good enough, stop
            for r = rhigh:-1:rlow
              Xr = round.(X,digits = r)
              adjust = n/sum(tril(Xr))
              Xr = adjust*Xr  # make sure sum of entries equals n
              objr = LR_obj(A,Xr)
              xWxr = xWnorm(W,Xr)
              PrimalRound = xWxr/(2*gam) + objr
              roundTri = FullTriangleCheck(Xr)
              roundedGap = (PrimalRound-DualQP)/DualQP
              R = xWxr/(2*gam*objr)
              if roundTri < ConTol && abs(roundedGap) < GapTol
                      roundR = r
                      roundconverge = true
                      break
              end
            end
            specialPrint = true
          end
        end

        # Save current progress
        push!(primals,PrimalQP)
        push!(duals,DualQP)
        push!(gaps,gap)
        push!(ConViolation,lastConCheck)
        push!(LPobjs,objective)

        # Round things to print out results nicely
        gapround = round(gap,digits = 5)
        PriRound = round(PrimalQP,digits = 5)
        DuRound = round(DualQP,digits = 5)
        ob = round(objective,digits = 3)
        time = round(triTime,digits = 3)
        tr = round(lastConCheck,digits = 5)
        Bty =  round(-BtDual/(gam),digits = 4)
        println("Iter $iter: Dual = $DuRound, Primal = $PriRound, gap = $gapround, ConVio = $tr, 3Loop = $time, Obj: $ob")
        open(filename,"a") do f
            write(f, "Iter $iter: Dual = $DuRound, Primal = $PriRound, gap = $gapround, ConVio = $tr, 3Loop = $time, Obj: $ob \n")
        end

        # Print something extra if you perform the entrywise rounding step
        if specialPrint
                PR = round(PrimalRound,digits = 5); rgap = round(roundedGap,digits = 5); rTri = round(roundTri,digits = 5)
                println("\t Rounded: \tNewPrimal = $PR, gap = $rgap, ConVio = $rTri")
                open(filename,"a") do f
                    write(f, "\t Rounded: \tNewPrimal = $PR, gap = $rgap, ConVio = $rTri \n")
                end
        end

        return tricheck, lastConCheck, objective, LPobjs, primals, duals, gaps, ConViolation, gap, roundconverge, roundR,R, roundedGap, roundTri,Bty
end

# Enforce the constraint \sum_ij x_ij = n, part of the Leighton Rao sparsest
# cut relaxation
function the_sum_constraint!(X::Matrix{Float64},W::Matrix{Float64},SumCorrection::Float64,InvSumW::Float64)

    n = size(X,1)
    # Correction step
    sumX = 0
    for i = 1:n-1
        for j = i+1:n
            X[j,i] = X[j,i] + SumCorrection/W[j,i]
            sumX += X[j,i]
        end
    end

    constant = (sumX - n)/InvSumW
    for i = 1:n-1
        for j = i+1:n
            X[j,i] = X[j,i] - constant/W[j,i]
        end
    end
    return constant

end

# Enforce constraint X_ij >= 0 for the sparsest cut relaxation
# P is the set of correction variables
# W is the weights matrix
function box_constraints!(X::Matrix{Float64},W::Matrix{Float64},P::Matrix{Float64})

    n = size(X,1)
    for i = 1:n-1
        for j = i+1:n
            X[j,i] -= P[j,i]/W[j,i]
            thetaIplus = -X[j,i]*W[j,i]
            if thetaIplus > 0
                X[j,i] = 0.0
                P[j,i] = thetaIplus
            else
                P[j,i] = 0.0
            end
        end
    end
end


# Enforce triangle inequality constraints
function cyclic_triple_loopX!(X::Matrix{Float64},W::Matrix{Float64},
  now_corrections::Vector{Tuple{Int64,Float64}},
    next_corrections::Vector{Tuple{Int64,Float64}})

  n = size(X,1)
  epsi = 0

  correctionsLength = length(now_corrections)
  nowInd = 1
  nextKey = now_corrections[nowInd]

  @inbounds for i = 1:n-2
    for  j = i+1:n-1
     Wij = W[j,i]

    for  k = j+1:n
    Xik = X[k,i]
    Xjk = X[k,j]
    Xij = X[j,i]
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

      X[j,i] = Xij + cor*(Wik*Wjk/denom)
      X[k,i] = Xik - cor*(Wij*Wjk/denom)
      X[k,j] = Xjk - cor*(Wik*Wij/denom)
      Xij = X[j,i]
      Xik = X[k,i]
      Xjk = X[k,j]
      # Move along in the list of triplets with corrections
      if nowInd < correctionsLength
              nowInd +=1
              nextKey = now_corrections[nowInd]
      end
    end

    mu = (Xij - Xjk - Xik)

    if mu > epsi
      denom = Wij*Wjk + Wik*Wij + Wjk*Wik

      X[j,i] = Xij - mu*(Wik*Wjk/denom)
      X[k,i] = Xik + mu*(Wij*Wjk/denom)
      X[k,j] = Xjk + mu*(Wik*Wij/denom)
      # Next time we see this triple we have to correct
      push!(next_corrections,(ijkKey,mu))
    end

    ### Done checking triangle i,j,k

    ### Check triangle i,k,j
    Xij = X[j,i]
    Xik = X[k,i]
    Xjk = X[k,j]

    ijkKey = (i-1)*n^2+(k-1)*n+j
    # First see if this is the next triangle with a nonzero correction variable
    if ijkKey == nextKey[1]

            cor = nextKey[2]

      denom = Wij*Wjk + Wik*Wij + Wjk*Wik
      X[j,i] = X[j,i] - cor*(Wik*Wjk/denom)
      X[k,i] = X[k,i] + cor*(Wij*Wjk/denom)
      X[k,j] = X[k,j] - cor*(Wik*Wij/denom)
      Xij = X[j,i]
      Xik = X[k,i]
      Xjk = X[k,j]
      # Move along in the list of triplets with corrections
      if nowInd < correctionsLength
              nowInd +=1
              nextKey = now_corrections[nowInd]
      end
    else
      Xij = X[j,i]
      Xik = X[k,i]
      Xjk = X[k,j]
    end
    mu = (-Xij - Xjk + Xik)
    if mu > epsi

      denom = Wij*Wjk + Wik*Wij + Wjk*Wik

      X[j,i] = Xij + mu*(Wik*Wjk/denom)
      X[k,i] = Xik - mu*(Wij*Wjk/denom)
      X[k,j] = Xjk + mu*(Wik*Wij/denom)

      # Next time we see this triple we have to correct
      push!(next_corrections,(ijkKey,mu))
    end
    ### Done checking triangle i,k,j

    ### Triangle j,k,i
    Xij = X[j,i]
    Xik = X[k,i]
    Xjk = X[k,j]
    ijkKey = (j-1)*n^2+(k-1)*n+i
    # First see if this is the next triangle with a nonzero correction variable
    if ijkKey == nextKey[1]

            cor = nextKey[2]

      denom = Wij*Wjk + Wik*Wij + Wjk*Wik

      X[j,i] = X[j,i] - cor*(Wik*Wjk/denom)
      X[k,i] = X[k,i] - cor*(Wij*Wjk/denom)
      X[k,j] = X[k,j] + cor*(Wik*Wij/denom)
      Xij = X[j,i]
      Xik = X[k,i]
      Xjk = X[k,j]
      # Move along in the list of triplets with corrections
      if nowInd < correctionsLength
              nowInd +=1
              nextKey = now_corrections[nowInd]
      end
    else
      Xij = X[j,i]
      Xik = X[k,i]
      Xjk = X[k,j]
    end

    mu = (-Xij + Xjk - Xik)

    if mu > epsi
      denom = Wij*Wjk + Wik*Wij + Wjk*Wik

      X[j,i] = Xij + mu*(Wik*Wjk/denom)
      X[k,i] = Xik + mu*(Wij*Wjk/denom)
      X[k,j] = Xjk - mu*(Wik*Wij/denom)

      # Next time we see this triple we have to correct
      push!(next_corrections,(ijkKey,mu))
    end
    ### Done checking triangle j,k,i

end
end
end # end triple for loop

end # end TripleLoop! function


# """
# LamCC_DandW:
# Builds the weights matrix W and dissimilarity matrix D
# for a LambdaCC LP relaxation problem
# """
function LamCC_DandW(A::SparseMatrixCSC{Float64,Int64},lam::Float64)
  n = size(A,1)
  D = zeros(Float64,n,n)
  W = (1-lam)*ones(n,n)
  for i = 1:n-1
    for j = i+1:n
      if A[i,j] < .1
        D[j,i] = 1
        W[j,i] = lam
      end
    end
  end
  D, W
end


# Builds the weights matrix and desired distance matrix D
# for a LambdaCC LP relaxation problem, this time degree-weighted
function LamCC_degweighted(A::SparseMatrixCSC{Float64,Int64},lam::Float64)
  n = size(A,1)
  D = zeros(Float64,n,n)
  W = zeros(n,n)
  d = sum(A,dims =2)
  for i = 1:n-1
      for j = i+1:n
          if A[i,j] < .1
              D[j,i] = 1
              W[j,i] = d[i]*d[j]*lam
          else
              if 1 - d[i]*d[j]*lam < 0
                  W[j,i] = d[i]*d[j]*lam - 1
                  D[j,i] = 1
              else
                  W[j,i] = (1-d[i]*d[j]*lam)
              end
          end
      end
  end
  return D, W
end

# """
# LeightonRaoQP_Initialize:
# Initialize solution vector and weights matrix for the Leighton Rao quadratic
# programming relaxation.
# """
function LeightonRaoQP_Initialize(A::SparseMatrixCSC{Float64,Int64},lam::Float64,gam::Float64)

  n = size(A,1)
  # X is a dense matrix
  X = zeros(Float64,n,n)
  W = lam*ones(n,n)
  @inbounds for i = 1:n-1
    for j = i+1:n
      if A[i,j] == 1
        X[j,i] = -gam
        W[j,i] = 1
      end
    end
  end
  for i = 1:n
    W[i,i] = 0.0
  end
  X, tril(W)
end


# """
# TriangleCheck:
# Returns whether or not all triangle inequality constraints are satisfied
# to within the desired tolerance. Returns as soon as it finds a constraint
# violation.
# """
function TriangleCheck(D::Matrix{Float64},tol::Float64)
  n = size(D,1)
  @inbounds for i = 1:n-2
    for j = i+1:n-1
      for k = j+1:n
        a = D[j,i]
        b = D[k,i]
        c = D[k,j]
        if a - b - c > tol || b - c - a > tol || c - a - b > tol
          return false
        end
      end
    end
  end
  return true
end

# FullTriangleCheck
# Returns the worst triangle violation in the whole matrix
function FullTriangleCheck(D::Matrix{Float64})
  n = size(D,1)
  maxi = 0.0
  @inbounds for i = 1:n-2
    for j = i+1:n-1
      a = D[j,i]
      for k = j+1:n
        b = D[k,i]
        c = D[k,j]

        #vio = maximum([a-b-c,b-a-c,c-a-b])
        if a-b > maxi && a-c > maxi && a-b-c > maxi
          maxi = a-b-c
        end
        if b - a > maxi && b-c > maxi && b-c-a > maxi
          maxi = b-c-a
        end
        if c-a > maxi && c-c > maxi && c-a-b > maxi
          maxi = c-a-b
        end

      end
    end
  end
  maxi
end

# """
# DoubleCheck:
# The L1 metric nearness problem (which the correlation clustering LP relaxation
# is a special case of) contains constraints eij - fij <= 0, -eij - fij <= 0.
# This function checks if they are satisfied to within a certain tolerance.
# """
function DoubleCheck(E::Matrix{Float64},F::Matrix{Float64},tol::Float64)
n = size(E,1)
for i = 1:n-1
  for j = i+1:n
    eij = E[j,i]
    fij = F[j,i]
    if eij - fij > tol || -eij - fij > tol
      return false
    end
  end
end
return true
end


# FullTriangleCheck
# Returns the worst triangle violation in the whole matrix
function FullTriangleCheck(D::Matrix{Float64})
  n = size(D,1)
  maxi = 0.0
  @inbounds for i = 1:n-2
    for j = i+1:n-1
      a = D[j,i]
      for k = j+1:n
        b = D[k,i]
        c = D[k,j]

        #vio = maximum([a-b-c,b-a-c,c-a-b])
        if a-b > maxi && a-c > maxi && a-b-c > maxi
          maxi = a-b-c
        end
        if b - a > maxi && b-c > maxi && b-c-a > maxi
          maxi = b-c-a
        end
        if c-a > maxi && c-c > maxi && c-a-b > maxi
          maxi = c-a-b
        end

      end
    end
  end
  maxi
end

# Checking all constraints for the correlation clustering LP problem
function FullConstraintCheck(X::Matrix{Float64},E::Matrix{Float64},F::Matrix{Float64})

  tri = FullTriangleCheck(X)
  n = size(X,1)
  maxi = 0.0

  for i = 1:n-1
    for j = i+1:n
      eij = E[j,i]
      fij = F[j,i]
      if eij - fij > maxi
        maxi = eij - fij
      end
      if -eij - fij > maxi
        maxi = -eij - fij
      end
    end
  end
  doublecheck = maxi

  return round(max(tri,doublecheck),digits = 4), round(tri,digits = 4)
end


# Evaluate the LambdaCC LP relaxed objective, given a distance matrix D
function LPcc_obj(A::SparseMatrixCSC{Float64,Int64},D::Matrix{Float64},lam::Float64)
  n = size(A,1)
  # assert(issymmetric(D))
  # assert(issymmetric(A))
  numedges = countnz(A)/2
  lccBound = sum((A[j,i]-lam)*D[j,i] for i=1:n-1 for j = i+1:n) + lam*(n*(n-1)/2 - numedges)
  return lccBound
end

# This is effectively a vector dot product, but for matrices.
# Specifically, this corresponds to the linear program objective score for
# variables F.
function LPobj(F::Matrix{Float64},W::Matrix{Float64})
  n = size(F,1)
  obj = 0.0
  for i = 1:n-1
    for j = i+1:n
      obj += F[j,i]*W[j,i]
    end
  end
  return obj
end

# Computes the norm for the vector of variables for the correlation clustering problem
function Wnorm(W::Matrix{Float64},E::Matrix{Float64},F::Matrix{Float64})

n = size(W,1)
out = 0.0
for i = 1:n-1
  for j = i+1:n
    out += (E[j,i]^2 + F[j,i]^2)*W[j,i]
  end
end
return out

end

# Computes the norm for the vector of variables for the Leighton-Rao relaxaton LP
function xWnorm(W::Matrix{Float64},X::Matrix{Float64})
  n = size(W,1)
  out = 0.0
  for i = 1:n-1
    for j = i+1:n
      out += (X[j,i]^2)*W[j,i]
    end
  end
  return out
end

# Computes the objective for the Leighton Rao sparsest cut LP relaxation
function LR_obj(A::SparseMatrixCSC{Float64,Int64},X::Matrix{Float64})
  n = size(A,1)
  lr = sum((A[j,i]*X[j,i] for i=1:n-1 for j = i+1:n))
  return lr
end
