#
# I've used several versions of these algorithms, and now I want to implement
# the cleanest and fastest version of Dykstra's possible to compare all future
# attempts against.
#
#   Coded by Nate Veldt on February 19, 2018

include("Tricon_Helper.jl")

import Base.hash
hash(x::Integer) = UInt64(x)

# DYKSTRA_LAMCC_TFA
#
# A Triangle Fixing Algorithm for the LambdaCC LP relaxation
#     based on Dykstra's projection algorithm
#
function Dykstra_lamCC_TFA(A::SparseMatrixCSC{Float64,Int64},tol::Float64=1e-3,TriTol::Float64=1e-3,
                lam::Float64=0.5,filename::String="DykstraLamCCoutput",gam::Float64=10.0,maxits::Int64=1000)

        n = size(A,1)
        open(filename, "w") do f
                write(f, "Output from Hildreths lamCC TFA\n")
                write(f, "Lambda = $lam, gamma = $gam, tol = $tol, TriTol = $TriTol \n")
        end

        E = zeros(n,n)
        F = -gam*ones(n,n)
        D,W = ConstructDandW(A)
        # # Matrices for double-loop correction terms (purposely dense)
        P = zeros(n,n)
        Q = zeros(n,n)

        # Correction term vector for triangle constraints
        current_corrections = Vector{Tuple{Int64,Float64}}()
        next_corrections = Vector{Tuple{Int64,Float64}}()
        # First time through constraints
        cyclic_triangle_constraints!(A,D,E,W,current_corrections,next_corrections)
        cyclic_double_constraints!(E,F,P,Q)

        iter = 0
        lastTri = 1.0
        objective = 0.0
        while iter < maxits
                iter += 1
                vecStart = [vec(E);vec(F)]

                # An empty vector to fill future correction terms in
                current_corrections = next_corrections
                next_corrections = Vector{Tuple{Int64,Float64}}()

                tic()
                cyclic_triangle_constraints!(A,D,E,W,current_corrections,next_corrections)
                TriTime = toq()

                cyclic_double_constraints!(E,F,P,Q)
                vecNext = [vec(E);vec(F)]
                vecChange = sqrt(mapreduce(z -> abs(z[1]-z[2])^2, +, zip(vecStart,vecNext)))
                #vecChange = norm(xStart-xNext)

                tricheck, lastTri, objective = report_progress(A,E+D,current_corrections,filename,vecChange,iter,TriTol,lam,lastTri,TriTime)

                if vecChange < Etol && tricheck
                  break
                end
        end
        X = E+D
        Finaltr = FullTriangleCheck(X)
        Finalobj = LPcc_obj(A,E+D,lam)
        Finalits = iter
        return X, Finaltr, Finalobj, Finalits
end

function report_progress(A::SparseMatrixCSC{Float64,Int64},X::Matrix{Float64},current_corrections::Vector{Tuple{Int64,Float64}},
    filename::String,xChange::Float64,iter::Int64,TriTol::Float64,lam::Float64,lastTri::Float64,triTime::Float64)

        tricheck = TriangleCheck(X,TriTol)
        objective = LPcc_obj(A,X,lam)
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


function cyclic_triangle_constraints!(A::SparseMatrixCSC{Float64,Int64},D::Matrix{Int8},E::Matrix{Float64},
    W::Matrix{Float64},current_corrections::Vector{Tuple{Int64,Float64}},next_corrections::Vector{Tuple{Int64,Float64}})

curLength = length(current_corrections)
curKey = 0
curInd = 0
if curLength > 0
        curInd = 1
        curKey = current_corrections[curInd][1]
end
n = size(A,1)
@inbounds for i = 1:n-2
    for j = i+1:n-1
        Wij = W[j,i]
        Dij = D[j,i]
        for k = j+1:n
            Wik = W[k,i]
            Wjk = W[k,j]
            Dik = D[k,i]
            Djk = D[k,j]

            Dijk = Float64(Dij - Dik - Djk)
            Djki = Float64(Djk - Dik - Dij)
            Dikj = Float64(Dik - Djk - Dij)

            denom = Wik*Wjk + Wij*Wik + Wij*Wjk
            # i,j,k here satisfies i < j < k
            # There are three ways to order these, each with a different constraint
            # curInd, curKey = ABCproj!(E,current_corrections,next_corrections,curInd,curLength,curKey,i,j,Wij,Wik,Wjk,n,denom,(i-1)*n^2+(j-1)*n + k,k,i,k,j,Dijk)
            # curInd, curKey = ABCproj!(E,current_corrections,next_corrections,curInd,curLength,curKey,i,k,Wik,Wij,Wjk,n,denom,(i-1)*n^2+(k-1)*n + j,j,i,k,j,Dikj)
            # curInd, curKey = ABCproj!(E,current_corrections,next_corrections,curInd,curLength,curKey,j,k,Wjk,Wij,Wik,n,denom,(j-1)*n^2+(k-1)*n + i,j,i,k,i,Djki)
            curInd, curKey = ABCproj!(E,current_corrections,next_corrections,curInd,curLength,curKey,i,j,Wij,Wik,Wjk,n,denom,(i-1)*n^2+(j-1)*n + k,k,i,k,j,Dijk)
            curInd, curKey = ABCproj!(E,current_corrections,next_corrections,curInd,curLength,curKey,i,k,Wik,Wij,Wjk,n,denom,(i-1)*n^2+(k-1)*n + j,j,i,k,j,Dikj)
            curInd, curKey = ABCproj!(E,current_corrections,next_corrections,curInd,curLength,curKey,j,k,Wjk,Wij,Wik,n,denom,(j-1)*n^2+(k-1)*n + i,j,i,k,i,Djki)

        end
    end
end

end

# This is the projection where we subtract the correction term and then project.
# It's how you most often see Dykstra's method presented, rather than the
# typical presentation of Hildreth's, even though they are equivalent.
function ABCproj!(E::Matrix{Float64},current_corrections::Vector{Tuple{Int64,Float64}},next_corrections::Vector{Tuple{Int64,Float64}},
    curInd::Int64,curLength::Int64,curKey::Int64,a::Int64,b::Int64,wab::Float64,wac::Float64,wbc::Float64,n::Int64,denom::Float64,
    ijkKey::Int64,maxAC::Int64,minAC::Int64,maxBC::Int64,minBC::Int64,Dabc::Float64)

    Eab = E[b,a]
    Eac = E[maxAC,minAC]
    Ebc = E[maxBC,minBC]

    # Fix: instead of updating curInd alone, also update Next Corr index.
    if ijkKey == curKey
        corr = current_corrections[curInd][2]
        if curInd < curLength
                curInd += 1
                curKey = current_corrections[curInd][1]
        end

        E[b,a] = Eab + corr*wbc*wac/denom
        E[maxAC,minAC] = Eac - corr*wbc*wab/denom
        E[maxBC,minBC] = Ebc - corr*wac*wab/denom

        Eab = E[b,a]
        Eac = E[maxAC,minAC]
        Ebc = E[maxBC,minBC]
    end

    delta = Eab - Eac - Ebc + Dabc
    if delta > 0
        E[b,a] = Eab - delta*wbc*wac/denom
        E[maxAC,minAC] = Eac + delta*wab*wbc/denom
        E[maxBC,minBC] = Ebc + delta*wab*wac/denom
        push!(next_corrections,(ijkKey,delta))
    end

    return curInd, curKey
end


function cyclic_double_constraints!(E::Matrix{Float64},F::Matrix{Float64},
    P::Matrix{Float64},Q::Matrix{Float64})

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
      P[j,i] = 0.0
    end

    #  E - F <= 0
    # corrections
    cor = Q[j,i]
    Eij = E[j,i] - cor
    Fij = F[j,i] + cor
    delta = Eij - Fij
    if delta > 0.0
      E[j,i] = Eij - delta/2
      F[j,i] = Fij + delta/2
      Q[j,i] = -delta/2
    else
      Q[j,i] = 0.0
    end

  end
end

end
