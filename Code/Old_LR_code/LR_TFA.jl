#
# I've used several versions of these algorithms, and now I want to implement
# the cleanest and fastest version of Dykstra's possible to compare all future
# attempts against.
#
#   Coded by Nate Veldt on February 19, 2018

include("Tricon_Helper.jl")

using JuMP
using Gurobi

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
        cVec = zeros(n,n)
        P = zeros(n,n)
        SumCorrection = 0.0
        InvSumWeights = 0.0
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
        combine_constraints!(X,cVec)
        @show sum(X), minimum(X), maximum(X)
        SumCorrection = cyclic_triple_loop!(X,W,current_corrections,next_corrections,SumCorrection,InvSumWeights)

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
                # xx = combine_constraints!(X,cVec)
                # X = tril(xx)
                # #@show sum(X), minimum(X), maximum(X)
                SumCorrection = cyclic_triple_loop!(X,W,current_corrections,next_corrections,SumCorrection,InvSumWeights)
                TriTime = toq()
                # xx = combine_constraints!(X,cVec)
                # X = tril(xx)

                vecNext = vec(X)
                vecChange = sqrt(mapreduce(z -> abs(z[1]-z[2])^2, +, zip(vecStart,vecNext)))
                #vecChange = norm(xStart-xNext)

                tricheck, lastTri, objective = report_progress(A,X,current_corrections,filename,vecChange,iter,TriTol,lam,lastTri,TriTime)

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

function combine_constraints!(X::Matrix{Float64},cVec::Matrix{Float64})

    X = X - cVec
    n = size(X,1)
    m = Model(solver=GurobiSolver(OutputFlag = 0))
    @variable(m, x[1:n,1:n])

    # Minimize the weight of disagreements (relaxation)

    @objective(m, Min, sum( (X[j,i] - x[j,i])^2 for i=1:n-1 for j = i+1:n ) )

    for i = 1:n-1
        for j = i+1:n
            @constraint(m,x[j,i] <= 1)
            @constraint(m,x[j,i] >= 0)
        end
    end

    for i = 1:n
        @constraint(m,x[i,i] == 0)
        for j = 1:n
            @constraint(m,x[i,j] == x[j,i])
        end
    end
    # for j = 1:n
    #     for i = j:n
    #             @constraint(m,x[i,j] == 0)
    #     end
    # end
    # Constraint: \sum_{ij} x_{ij} = k for some value k, we'll use 1
    # but in Luca Trevisan's notes he uses something else because the definition
    # of sparsest cut he is working with is a multiplicative factor different

    @constraint(m,sum(x[j,i] for i=1:n-1 for j = i+1:n) == n)
    solve(m)
    xx = getvalue(x)
    return xx
    # @show sum(tril(xx)), minimum(tril(xx)), maximum(tril(xx))
    # X = tril(xx)
    cVec = xx - X

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


function cyclic_triple_loop!(X::Matrix{Float64},W::Matrix{Float64},now_corrections::Vector{Tuple{Int64,Float64}},next_corrections::Vector{Tuple{Int64,Float64}},SumCorrection::Float64,InvSumW::Float64)
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
#SumCorrection = sum_constraint!(X,W,SumCorrection,InvSumW)

return SumCorrection
end # end TripleLoop! function
