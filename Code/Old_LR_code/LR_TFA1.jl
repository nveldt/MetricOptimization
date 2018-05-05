#
# I've used several versions of these algorithms, and now I want to implement
# the cleanest and fastest version of Dykstra's possible to compare all future
# attempts against.
#
#   Coded by Nate Veldt on February 19, 2018

include("Tricon_Helper.jl")

import Base.hash
hash(x::Integer) = UInt64(x)

function ConstructDandW(A::SparseMatrixCSC{Float64,Int64},lam::Float64)
 n = size(A,1)
  Amat = zeros(Float64,n,n)
  W = lam*ones(n,n)
  for i = 1:n-1
      for j = i+1:n
          if A[i,j] > .1
              Amat[j,i] = 1
              W[j,i] = 1
          end
      end
  end
  return Amat, W
end

# DYKSTRA_LeightonRao_TFA
#
# A Triangle Fixing Algorithm for the Leighton-Rao Relaxation of sparsest cut
#
function Dykstra_LeightonRao_TFA(A::SparseMatrixCSC{Float64,Int64},Xtol::Float64=1e-3,TriTol::Float64=1e-5,
                lam::Float64=0.1,filename::String="DykstraLeightonRaoOutput",gam::Float64=100.0,maxits::Int64=1000)

        n = size(A,1)
        open(filename, "w") do f
                write(f, "Output from Hildreths lamCC TFA\n")
                write(f, "Lambda = $lam, gamma = $gam, tol = $Xtol, TriTol = $TriTol \n")
        end

        Amat,W = ConstructDandW(A,lam)
        X = -gamma*Amat
        #W = ones(n,n)
        # # Matrices for double-loop correction terms (purposely dense)
        P = zeros(n,n)
        SumCorrection = 0.0
        InvSumWeights = 0.0
        SlackMax = 0.0
        TotalSlack = 0.0
        for i = 1:n-1
            for j = i+1:n
                InvSumWeights += 1/W[j,i]
            end
        end
        @show InvSumWeights, n*(n-1)/2
        # Correction term vector for triangle constraints
        current_corrections = Vector{Tuple{Int64,Float64}}()
        next_corrections = Vector{Tuple{Int64,Float64}}()

        push!(current_corrections,(0,0.0))

        # First time through constraints
        box_constraints!(X,W,P)
        SumCorrection = sum_constraint!(X,W,SumCorrection,InvSumWeights)
        cyclic_triple_loop!(X,W,current_corrections,next_corrections)

        iter = 0
        lastTri = 1.0
        objective = 0.0
        while iter < maxits
                iter += 1
                vecStart = vec(X)

                # An empty vector to fill future correction terms in
                current_corrections = next_corrections
                next_corrections = Vector{Tuple{Int64,Float64}}()
                push!(current_corrections,(0,0.0))

                tic()
                box_constraints!(X,W,P)
                SumCorrection = sum_constraint!(X,W,SumCorrection,InvSumWeights)
                cyclic_triple_loop!(X,W,current_corrections,next_corrections)
                TriTime = toq()

                vecNext = vec(X)
                vecChange = sqrt(mapreduce(z -> abs(z[1]-z[2])^2, +, zip(vecStart,vecNext)))
                #vecChange = norm(xStart-xNext)

                tricheck, lastTri, objective, SlackMax, TotalSlack, SumCorrection = new_report_progress(A,X,P,current_corrections,filename,vecChange,iter,TriTol,lam,lastTri,SlackMax,TotalSlack,TriTime,SumCorrection)

                if vecChange < Xtol && tricheck && iter > 100
                  break
                end
        end
        Finaltr = FullTriangleCheck(X)
        Finalobj = LR_obj(A,X)
        Finalits = iter
        return X, Finaltr, Finalobj, Finalits
end

function report_progress(A::SparseMatrixCSC{Float64,Int64},X::Matrix{Float64},current_corrections::Vector{Tuple{Int64,Float64}},
    filename::String,xChange::Float64,iter::Int64,TriTol::Float64,lam::Float64,lastTri::Float64,triTime::Float64)

        tricheck = TriangleCheck(X,TriTol)
        objective = LR_obj(A,X)
        ch = round(xChange,3)
        ob = round(objective,3)
        time = round(triTime,3)
        tr = round(lastTri,3)
        if iter%20 == 1
          lastTri = FullTriangleCheck(X)
        end
        nnzDelta = length(current_corrections)
        println("Iteration $iter, ||x_k - x_{k+1}|| = $ch, 3Loop = $time, Obj: $ob, Last TriCheck = $tr, nnz(delta) = $nnzDelta")

        open(filename, "a") do f
          write(f,"Iteration $iter, ||x_k - x_{k+1}|| = $ch, 3Loop = $time, Obj: $ob, Last TriCheck = $tr\n")
        end

        return tricheck, lastTri, objective
end

# Check the worst violation in the KKT conditions
function CheckKKT(X::Matrix{Float64},now_corrections::Vector{Tuple{Int64,Float64}},P::Matrix{Float64},SumCorrection::Float64)

    n = size(X,1)
    correctionsLength = length(now_corrections)
    nowInd = 1
    # Grab next triplet in the list
    nextKey = now_corrections[nowInd]

    maxi = 0.0
    SlackMax = 0.0
    TotalSlack = 0.0

    @inbounds for i = 1:n-2
      for j = i+1:n-1
        a = X[j,i]
        for k = j+1:n
          b = X[k,i]
          c = X[k,j]

          dijk = a - b - c
          dikj = b - c - a
          djki = c - a - b

          # Check for complementary slackness

          ## Check triangle i,j,k
         ijkKey = (i-1)*n^2+(j-1)*n+k
         if ijkKey == nextKey[1]

            cor = nextKey[2]
            slackness = abs(cor*dijk)
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
        # Now check for the slackness condition on the sum constraint
        sumX = sum(X)
        violation = (sumX - n)                      #should be zero
        slackness = abs(violation*SumCorrection)
        TotalSlack += slackness
        if slackness > SlackMax
            SlackMax = slackness
        end

        # Now check whether x >= 0
        for i = 1:n-1
            for j = i+1:n
                violation = X[j,i]
                slackness = abs(violation*P[j,i])
                TotalSlack += slackness
                if slackness > SlackMax
                    SlackMax = slackness
                end
            end
        end

    return round(maxi,3), round(SlackMax,3), round(TotalSlack,3)
  end

function new_report_progress(A::SparseMatrixCSC{Float64,Int64},X::Matrix{Float64},P::Matrix{Float64},current_corrections::Vector{Tuple{Int64,Float64}},
    filename::String,xChange::Float64,iter::Int64,TriTol::Float64,lam::Float64,lastTri::Float64,SlackMax::Float64,TotalSlack::Float64,triTime::Float64,SumCorrection::Float64)

        tricheck = TriangleCheck(X,TriTol)
        objective = LR_obj(A,X)
        ch = round(xChange,3)
        ob = round(objective,3)
        time = round(triTime,3)
        tr = round(lastTri,3)
        if iter%10 == 1
          @time lastTri = FullTriangleCheck(X)
          @time lastTri, SlackMax, TotalSlack = CheckKKT(X,current_corrections,P,SumCorrection)
        end
        nnzDelta = length(current_corrections)
        println("Iteration $iter, ||x_k - x_{k+1}|| = $ch, 3Loop = $time, Obj: $ob, ConVio = $lastTri, SlackVio = $SlackMax , T-Slack = $TotalSlack")

        open(filename, "a") do f
          write(f,"Iteration $iter, ||x_k - x_{k+1}|| = $ch, 3Loop = $time, Obj: $ob, Last TriCheck = $tr\n")
        end

        return tricheck, lastTri, objective, SlackMax, TotalSlack, SumCorrection
end

function sum_constraint!(X::Matrix{Float64},W::Matrix{Float64},SumCorrection::Float64,InvSumW::Float64)

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

function box_constraints!(X::Matrix{Float64},W::Matrix{Float64},P::Matrix{Float64})

    n = size(X,1)
    for i = 1:n-1
        for j = i+1:n
            X[j,i] -= P[j,i]

            Xij = X[j,i]
            if Xij < 0
                X[j,i] = 0
                P[j,i] = -Xij
            else
                P[j,i] = 0.0
            end

        end
    end
end

function cyclic_triple_loop!(X::Matrix{Float64},W::Matrix{Float64},now_corrections::Vector{Tuple{Int64,Float64}},next_corrections::Vector{Tuple{Int64,Float64}})
  n = size(X,1)
  epsi = 0

  correctionsLength = length(now_corrections)
  nowInd = 1
  # Grab next triplet in the list
  nextKey = now_corrections[nowInd]

  ##nextTriplet = shift!(now_corrections)

  @inbounds for i = 1:n-2
    for  j = i+1:n-1
     Xij = X[j,i]
     Wij = W[j,i]

    for  k = j+1:n
    Xik = X[k,i]
    Xjk = X[k,j]
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
