#
# This is a newer implementation of Hildreth's algorithm (the special case of
# Dykstra's algorithm for half-space constraints)
#
# Major differences incldue:
#   1. We operate directly on distances X = (x_ij) rather than E = (e_ij)
#   2. This uses a dictionary TripletDelta to go from indices i,j,k to correction
#           variable delta (and we map i,j,k to a single unique index)
#   3. We include a parameter alpha that ranges between 0 and 2. Previous
#           previous implementations had alpha = 1 fixed.
#
#   Coded by Nate Veldt on February 19, 2018

include("Tricon_Helper.jl")

# HILDRETH_LAMCC_TFA
#
# A Triangle Fixing Algorithm for the LambdaCC LP relaxation
#     based on Hildreth's projection algorithm
#
# This solves the optimization problem:
#
# min sum_{i< j} w_ij y_ij + 1/(gamma) \sum_{i<j} w_ij (y_ij)^2 (double check that...)
#
# such that x_ij \leq x_ik + x_jk for all triplets i,j,k
#           x_ij - y_ij \leq d_ij
#           -x_ij - y_ij \leq -d_ij for all pairs i,j
#
function Hildreth_lamCC_TFA(A::SparseMatrixCSC{Float64,Int64},Xtol::Float64,TriTol::Float64,
                lam::Float64,filename::String,gam::Float64,maxits::Int64,alpha::Float64)

        n = size(A,1)
        open(filename, "w") do f
                write(f, "Output from Hildreths lamCC TFA\n")
                write(f, "Lambda = $lam, gamma = $gam, Xtol = $Xtol, TriTol = $TriTol, alpha = $alpha \n")
        end

        X = zeros(n,n)
        D = zeros(Int8,n,n)
        for i = 1:n-1
            for j = i+1:n
                if A[i,j] < .1
                    D[j,i] = 1
                    X[j,i] = 1.0
                end
            end
        end
        W = lam*ones(n,n)
        for i = 1:n-1
            for j = i+1:n
                if A[i,j] > .1
                    W[j,i] = 1-lam
                end
            end
        end

        #X = zeros(n,n)
        # Initial vector
        Y = -gam*ones(n,n)

        # Matrices for double-loop correction terms (purposely dense)
        P = zeros(n,n)
        Q = zeros(n,n)

        # Correction terms for triangle constraints
        triplet_corrections = Dict{Int64,Float64}()
        sizehint!(triplet_corrections,n^2)

        iter = 0
        lastTri = 1.0

        while iter < maxits
                iter += 1
                xStart = [vec(X);vec(Y)]

                tic()
                noncyclic_triangle_constraints!(A,D,X,W,triplet_corrections,alpha,Int64(ceil(n^2)))         #xij <= xik + xjk
                TriTime = toq()

                cyclic_double_constraints!(D,X,Y,P,Q,W,alpha)                           #xij - yij <= dij, etc.

                xNext = [vec(X);vec(Y)]
                xChange = norm(xStart-xNext)

                tricheck, lastTri = report_progress(A,X,triplet_corrections,filename,xChange,iter,TriTol,lam,lastTri,TriTime)

                if xChange < Etol && tricheck
                  break
                end
        end

        tic()
        for i = 1:n-2
            for j = i+1:n-2
                for k = j+1:n
                    ijkKey = (i-1)*n^2+(j-1)*n + k
                    corr = get(triplet_corrections,ijkKey,0.0)

                    # if corr == 0
                    #     delete!(triplet_corrections,ijkKey)
                    # else
                    #     triplet_corrections[ijkKey] = corr
                    # end
                end
            end
        end
        triloop = toc()

        @show triloop

        lastTri = FullTriangleCheck(X)
        return X, lastTri
end

function report_progress(A::SparseMatrixCSC{Float64,Int64},X::Matrix{Float64},triplet_corrections::Dict{Int64,Float64},
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
        nnzDelta = length(triplet_corrections)
        println("Iteration $iter, ||x_k - x_{k+1}|| = $ch, 3Loop = $time, Obj: $ob, Last TriCheck = $tr, nnz(delta) = $nnzDelta")

        open(filename, "a") do f
          write(f,"Iteration $iter, ||x_k - x_{k+1}|| = $ch, 3Loop = $time, Obj: $ob, Last TriCheck = $tr\n")
        end

        return tricheck, lastTri
end


function cyclic_triangle_constraints!(A::SparseMatrixCSC{Float64,Int64},D::Matrix{Int8},X::Matrix{Float64},
    W::Matrix{Float64},triplet_corrections::Dict{Int64,Float64},alpha::Float64)

n = size(A,1)
@inbounds for i = 1:n-2
    for j = i+1:n-1
        Wij = W[j,i]
        for k = j+1:n
            Wik = W[k,i]
            Wjk = W[k,j]
            denom = Wik*Wjk + Wij*Wik + Wij*Wjk
            # i,j,k here satisfies i < j < k
            # There are three ways to order these, each with a different constraint
            ABCproj!(X,triplet_corrections,alpha,i,j,Wij,Wik,Wjk,n,denom,(i-1)*n^2+(j-1)*n + k,k,i,k,j)
            ABCproj!(X,triplet_corrections,alpha,i,k,Wik,Wij,Wjk,n,denom,(i-1)*n^2+(k-1)*n + j,j,i,k,j)
            ABCproj!(X,triplet_corrections,alpha,j,k,Wjk,Wij,Wik,n,denom,(j-1)*n^2+(k-1)*n + i,j,i,k,i)
        end
    end
end

end

# This visits triangles at random and performs projections.
# It will converge, but it is slow (perhaps because of the time it takes to
# generate random numbers)
function noncyclic_triangle_constraints!(A::SparseMatrixCSC{Float64,Int64},D::Matrix{Int8},X::Matrix{Float64},
    W::Matrix{Float64},triplet_corrections::Dict{Int64,Float64},alpha::Float64,numtimes::Int64)

n = size(A,1)
RR = rand(1:n,numtimes,3)
Threads.@threads for t = 1:numtimes
    i,j,k = sort(RR[t,:])
    Wij = W[j,i]
    Wik = W[k,i]
    Wjk = W[k,j]
    denom = Wik*Wjk + Wij*Wik + Wij*Wjk
    # i,j,k here satisfies i < j < k
    # There are three ways to order these, each with a different constraint
    ABCproj!(X,triplet_corrections,alpha,i,j,Wij,Wik,Wjk,n,denom,(i-1)*n^2+(j-1)*n + k,k,i,k,j)
    ABCproj!(X,triplet_corrections,alpha,i,k,Wik,Wij,Wjk,n,denom,(i-1)*n^2+(k-1)*n + j,j,i,k,j)
    ABCproj!(X,triplet_corrections,alpha,j,k,Wjk,Wij,Wik,n,denom,(j-1)*n^2+(k-1)*n + i,j,i,k,i)
end

end


function ABCproj!(X::Matrix{Float64},triplet_corrections::Dict{Int64,Float64},alpha::Float64,
    a::Int64,b::Int64,wab::Float64,wac::Float64,wbc::Float64,n::Int64,denom::Float64,ijkKey::Int64,AC1::Int64,AC2::Int64,BC1::Int64,BC2::Int64)

    # Project onto the constraint convex set defined by xab <= xac + xbc
    # where it is NOT necessarily the case that a < b < c
    corr = pop!(triplet_corrections,ijkKey,0.0)
    Xab = X[b,a]
    Xac = X[AC1,AC2]
    Xbc = X[BC1,BC2]

    delta = Xab - Xac - Xbc

    noncorr = -alpha*delta*wab*wac*wbc/denom
    if corr < noncorr
        X[b,a] += corr/wab
        X[BC1,BC2] -= corr/wbc
        X[AC1,AC2] -= corr/wac
    else
        X[b,a] += noncorr/wab
        X[BC1,BC2] -= noncorr/wbc
        X[AC1,AC2] -= noncorr/wac
        triplet_corrections[ijkKey] = corr - noncorr
    end

end


# This is the projection where we subtract the correction term and then project.
# It's how you most often see Dykstra's method presented, rather than the
# typical presentation of Hildreth's, even though they are equivalent.
function ABCproj2!(X::Matrix{Float64},triplet_corrections::Dict{Int64,Float64},alpha::Float64,
    a::Int64,b::Int64,wab::Float64,wac::Float64,wbc::Float64,n::Int64,denom::Float64,ijkKey::Int64,AC1::Int64,AC2::Int64,BC1::Int64,BC2::Int64)

    corr = pop!(triplet_corrections,ijkKey,0.0)
    Xab = X[b,a]
    Xac = X[AC1,AC2]
    Xbc = X[BC1,BC2]
    # Correction step

    if corr > 0
        X[b,a] = Xab + corr*wbc*wac/denom
        X[AC1,AC2] = Xac - corr*wbc*wab/denom
        X[BC1,BC2] = Xbc - corr*wac*wab/denom

        Xab = X[b,a]
        Xac = X[AC1,AC2]
        Xbc = X[BC1,BC2]
    end

    delta = Xab - Xac - Xbc
    if delta > 0
        X[b,a] = Xab - alpha*delta*wbc*wac/denom
        X[AC1,AC2] = Xac + alpha*delta*wab*wbc/denom
        X[BC1,BC2] = Xbc + alpha*delta*wab*wac/denom
        triplet_corrections[ijkKey] = delta
    end

end


function cyclic_triangle_constraints2!(A::SparseMatrixCSC{Float64,Int64},D::Matrix{Int8},X::Matrix{Float64},
    W::Matrix{Float64},triplet_corrections::Dict{Int64,Float64},alpha::Float64)

n = size(A,1)
@inbounds for i = 1:n-2
    for j = i+1:n-1
        wij = W[j,i]
        for k = j+1:n
            wik = W[k,i]
            wjk = W[k,j]
            denom = wij*wjk + wik*wij + wjk*wik

            # i,j,k here satisfies i < j < k
            # There are three ways to order these, each with a different constraint

            ## Check Triangle i,j,k
            ijkKey = (i-1)*n^2+(j-1)*n + k
            corr = pop!(triplet_corrections,ijkKey,0.0)

            if corr > 0
                # Correction step
                X[j,i]  += corr*wik*wjk/denom
                X[k,i]  -= corr*wij*wjk/denom
                X[k,j]  -= corr*wij*wik/denom
            end

            # Check constraint
            delta = X[j,i] - X[k,i] - X[k,j]
            if delta > 1e-8
                X[j,i] -= delta*wik*wjk/denom
                X[k,i] += delta*wij*wjk/denom
                X[k,j] += delta*wij*wik/denom
                triplet_corrections[ijkKey] = delta
            end

            ## Check Triangle i,k,j
            ijkKey = (i-1)*n^2+(k-1)*n + j
            corr = pop!(triplet_corrections,ijkKey,0.0)

            if corr > 0
            # Correction step
            X[j,i]  -= corr*wik*wjk/denom
            X[k,i]  += corr*wij*wjk/denom
            X[k,j]  -= corr*wij*wik/denom
            end

            # Check constraint
            delta = -X[j,i] + X[k,i] - X[k,j]
            if delta > 1e-8
                X[j,i] += delta*wik*wjk/denom
                X[k,i] -= delta*wij*wjk/denom
                X[k,j] += delta*wij*wik/denom
                triplet_corrections[ijkKey] = delta
            end


            ## Check Triangle j,k,i
            ijkKey = (j-1)*n^2+(k-1)*n + i
            corr = pop!(triplet_corrections,ijkKey,0.0)

            if corr > 0
            # Correction step
            X[j,i]  -= corr*wik*wjk/denom
            X[k,i]  -= corr*wij*wjk/denom
            X[k,j]  += corr*wij*wik/denom
            end

            # Check constraint
            delta = -X[j,i] - X[k,i] + X[k,j]
            if delta > 1e-8
                X[j,i] += delta*wik*wjk/denom
                X[k,i] += delta*wij*wjk/denom
                X[k,j] -= delta*wij*wik/denom
                triplet_corrections[ijkKey] = delta
            end


        end
    end
end

end

function cyclic_double_constraints!(D::Matrix{Int8},X::Matrix{Float64},Y::Matrix{Float64},
    P::Matrix{Float64},Q::Matrix{Float64},W::Matrix{Float64},alpha::Float64)

    n = size(P,1)
    @inbounds for i = 1:n-1
      for j = i+1:n

        #  -X - Y <= -D
        cor = P[j,i]
        Xij = X[j,i] - cor
        Yij = Y[j,i] - cor
        delta = -Xij - Yij + D[j,i]

        if delta > 0
          X[j,i] = Xij + delta/2
          Y[j,i] = Yij + delta/2
          P[j,i] = delta/2
        else
          P[j,i] = 0.0
        end

        #  X - Y <= -D
        cor = Q[j,i]
        Xij = X[j,i] + cor
        Yij = Y[j,i] - cor
        delta = Xij - Yij - D[j,i]
        if delta > 0
          X[j,i] = Xij - delta/2
          Y[j,i] = Yij + delta/2
          Q[j,i] = delta/2
        else
          Q[j,i] = 0.0
        end

      end
    end

end

# This will converge to a stationary point, but isn't a minimum norm projection
function UncorrectedProjection(X::Matrix{Float64},triplet_corrections::Dict{Int64,Float64},alpha::Float64,
    a::Int64,b::Int64,wab::Float64,wac::Float64,wbc::Float64,n::Int64,denom::Float64,ijkKey::Int64,AC1::Int64,AC2::Int64,BC1::Int64,BC2::Int64)

    Xab = X[b,a]
    Xac = X[AC1,AC2]
    Xbc = X[BC1,BC2]

    delta = Xab - Xac - Xbc

    if delta > 0
        noncorr = delta*wab*wac*wbc/denom
        X[b,a] -= noncorr/wab
        X[BC1,BC2] += noncorr/wbc
        X[AC1,AC2] += noncorr/wac
        #triplet_corrections[ijkKey] = corr - noncorr
        return true
    else
        return false
    end
end
