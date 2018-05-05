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

        InvSumWeights = 0.0
        for i = 1:n-1
            for j = i+1:n
                InvSumWeights += 1/W[j,i]
            end
        end

        con = 1/gamma
        iter = 0
        lastTri = 1.0
        objective = 0.0
        while iter < maxits
                iter += 1
                vecStart = vec(X)

                tic()
                box_constraints!(X)
                cyclic_triple_loop!(X,W)
                sum_constraint!(X,W,InvSumWeights)
                TriTime = toq()


                Kcounter = iter
                sigmaK = con/(Kcounter+1)

                for i = 1:n-1
                  for j = i+1:n
                    X[j,i] = -gamma*sigmaK*Amat[j,i] + (1-sigmaK)*X[j,i]
                  end
                end

                vecNext = vec(X)
                vecChange = sqrt(mapreduce(z -> abs(z[1]-z[2])^2, +, zip(vecStart,vecNext)))
                #vecChange = norm(xStart-xNext)

                tricheck, lastTri, objective = report_progress(A,X,filename,vecChange,iter,TriTol,lam,lastTri,TriTime,gamma,Amat)

                if vecChange < Xtol && tricheck && iter > 100
                  break
                end
        end
        Finaltr = FullTriangleCheck(X)
        Finalobj = LR_obj(A,X)
        Finalits = iter
        return X, Finaltr, Finalobj, Finalits
end

function report_progress(A::SparseMatrixCSC{Float64,Int64},X::Matrix{Float64},
    filename::String,xChange::Float64,iter::Int64,TriTol::Float64,lam::Float64,lastTri::Float64,triTime::Float64,gamma::Float64,Amat::Matrix{Float64})

        tricheck = TriangleCheck(X,TriTol)
        objectiveLP = LR_obj(A,X)
        objective = norm(X + gamma*Amat)
        ch = round(xChange,4)
        ob = round(objective,4)
        oblp = round(objectiveLP,4)
        time = round(triTime,4)
        normX = round(norm(X),4)
        tr = round(lastTri,4)
        if iter%20 == 1
          lastTri = FullTriangleCheck(X)
        end
        println("Iteration $iter, ||x_k - x_{k+1}|| = $ch, 3Loop = $time, QPobj: $ob, LPobj: $oblp, normX: $normX Last TriCheck = $tr")

        open(filename, "a") do f
          write(f,"Iteration $iter, ||x_k - x_{k+1}|| = $ch, 3Loop = $time, Obj: $ob, Last TriCheck = $tr\n")
        end

        return tricheck, lastTri, objective
end

function sum_constraint!(X::Matrix{Float64},W::Matrix{Float64},InvSumW::Float64)
    n = size(X,1)
    # Correction step
    sumX = 0
    for i = 1:n-1
        for j = i+1:n
            sumX += X[j,i]
        end
    end
    constant = (sumX - n)/InvSumW
    for i = 1:n-1
        for j = i+1:n
            X[j,i] = X[j,i] - constant/W[j,i]
        end
    end
end

function box_constraints!(X::Matrix{Float64})
    n = size(X,1)
    for i = 1:n-1
        for j = i+1:n
            Xij = X[j,i]
            if Xij < 0
                X[j,i] = 0
            end
        end
    end
end


function cyclic_triple_loop!(X::Matrix{Float64},W::Matrix{Float64})
  n = size(X,1)

  ##nextTriplet = shift!(now_corrections)

  @inbounds for i = 1:n-2
    for  j = i+1:n-1
     Wij = W[j,i]
    for  k = j+1:n

    Xij = X[j,i]
    Xik = X[k,i]
    Xjk = X[k,j]
    Wik = W[k,i]
    Wjk = W[k,j]

    ### Check triangle i,j,k
    mu = (Xij - Xjk - Xik)
    if mu > 0
      denom = Wij*Wjk + Wik*Wij + Wjk*Wik
      X[j,i] = Xij - mu*(Wik*Wjk/denom)
      X[k,i] = Xik + mu*(Wij*Wjk/denom)
      X[k,j] = Xjk + mu*(Wik*Wij/denom)
    end

    ### Done checking triangle i,j,k

    ### Check triangle i,k,j
    mu = (-Xij - Xjk + Xik)
    if mu > 0

      denom = Wij*Wjk + Wik*Wij + Wjk*Wik
      X[j,i] = Xij + mu*(Wik*Wjk/denom)
      X[k,i] = Xik - mu*(Wij*Wjk/denom)
      X[k,j] = Xjk + mu*(Wik*Wij/denom)
    end
    ### Done checking triangle i,k,j

    mu = (-Xij + Xjk - Xik)
    if mu > 0
      denom = Wij*Wjk + Wik*Wij + Wjk*Wik
      X[j,i] = Xij + mu*(Wik*Wjk/denom)
      X[k,i] = Xik + mu*(Wij*Wjk/denom)
      X[k,j] = Xjk - mu*(Wik*Wij/denom)
    end
    ### Done checking triangle j,k,i

end
end
end # end triple for loop

end # end TripleLoop! function
