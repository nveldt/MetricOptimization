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
        current_corrections = Vector{Tuple{Int64,Int64,Int64,Float64}}()
        next_corrections = Vector{Tuple{Int64,Int64,Int64,Float64}}()

        push!(current_corrections,(0,0,0,0.0))
        # First time through constraints
        cyclic_triple_loop!(D,E,W,current_corrections,next_corrections)
        cyclic_double_loop!(E,F,P,Q)

        iter = 0
        lastTri = 1.0
        objective = 0.0
        while iter < maxits
                iter += 1
                vecStart = [vec(E);vec(F)]

                # An empty vector to fill future correction terms in
                current_corrections = next_corrections
                next_corrections = Vector{Tuple{Int64,Int64,Int64,Float64}}()

                #push!(current_corrections,(0,0,0,0.0))

                tic()
                cyclic_triple_loop!(D,E,W,current_corrections,next_corrections)
                TriTime = toq()

                cyclic_double_loop!(E,F,P,Q)
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

function report_progress(A::SparseMatrixCSC{Float64,Int64},X::Matrix{Float64},current_corrections::Vector{Tuple{Int64,Int64,Int64,Float64}},
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



function cyclic_double_loop!(E::Matrix{Float64},F::Matrix{Float64},
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

function cyclic_triple_loop!(D::Matrix{Int8},E::Matrix{Float64},W::Matrix{Float64},now_corrections::Vector{Tuple{Int64,Int64,Int64,Float64}},next_corrections::Vector{Tuple{Int64,Int64,Int64,Float64}})
  n = size(D,1)
  epsi = 1e-8

  correctionsLength = length(now_corrections)
  nowInd = 1;
  # Grab next triplet in the list
  nextTriplet = now_corrections[nowInd]

  ##nextTriplet = shift!(now_corrections)

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

    ### Check triangle i,j,k

    # First see if this is the next triangle with a nonzero correction variable
    if i == nextTriplet[1] && j == nextTriplet[2] && k == nextTriplet[3]
      cor = nextTriplet[4]

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
              nextTriplet = now_corrections[nowInd]
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
      push!(next_corrections,(i,j,k,mu))
    end

    ### Done checking triangle i,j,k


    ### Check triangle i,k,j
    if i == nextTriplet[1] && k == nextTriplet[2] && j == nextTriplet[3]
      cor = nextTriplet[4]

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
              nextTriplet = now_corrections[nowInd]
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
      push!(next_corrections,(i,k,j,mu))
    end
    ### Done checking triangle i,k,j

    ### Check triangle j,k,i
    if j == nextTriplet[1] && k == nextTriplet[2] && i == nextTriplet[3]
      cor = nextTriplet[4]

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
              nextTriplet = now_corrections[nowInd]
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
      push!(next_corrections,(j,k,i,mu))
    end
    ### Done checking triangle j,k,i

end
end
end # end triple for loop

end # end TripleLoop! function

## I have experimented with "cleaner" ways of laying out the triple loop projections
# However, it seems faster to just unroll the 3 constraints and write them all out
# explicitly.
function cyclic_triangle_constraints!(D::Matrix{Int8},E::Matrix{Float64},
    W::Matrix{Float64},current_corrections::Vector{Tuple{Int64,Int64,Int64,Float64}},next_corrections::Vector{Tuple{Int64,Int64,Int64,Float64}})

curLength = length(current_corrections)
curInd = 1
n = size(D,1)
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
            curInd = ABCproj!(E,current_corrections,next_corrections,curLength,curInd,i,j,k,Wij,Wik,Wjk,n,denom,k,i,k,j,Dijk)
            curInd = ABCproj!(E,current_corrections,next_corrections,curLength,curInd,i,k,j,Wik,Wij,Wjk,n,denom,j,i,k,j,Dikj)
            curInd = ABCproj!(E,current_corrections,next_corrections,curLength,curInd,j,k,i,Wjk,Wij,Wik,n,denom,j,i,k,i,Djki)
        end
    end
end

end

# This is the projection where we subtract the correction term and then project.
# It's how you most often see Dykstra's method presented, rather than the
# typical presentation of Hildreth's, even though they are equivalent.
function ABCproj!(E::Matrix{Float64},current_corrections::Vector{Tuple{Int64,Int64,Int64,Float64}},next_corrections::Vector{Tuple{Int64,Int64,Int64,Float64}},
    curLength::Int64,curInd::Int64,a::Int64,b::Int64,c::Int64,wab::Float64,wac::Float64,wbc::Float64,n::Int64,denom::Float64,
    maxAC::Int64,minAC::Int64,maxBC::Int64,minBC::Int64,Dabc::Float64)

    Eab = E[b,a]
    Eac = E[maxAC,minAC]
    Ebc = E[maxBC,minBC]

    # Fix: instead of updating curInd alone, also update Next Corr index.
    if a == current_corrections[curInd][1] && b == current_corrections[curInd][2] && c == current_corrections[curInd][3]
        corr = current_corrections[curInd][4]
        if curInd < curLength
                curInd += 1
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
        push!(next_corrections,(a,b,c,delta))
    end

    return curInd
end
