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
                cyclic_triangle_constraints!(A,D,X,Y,W,triplet_corrections,P,Q)         #xij <= xik + xjk
                TriTime = toc()

                #cyclic_double_constraints!(D,X,Y,P,Q,W,alpha)                           #xij - yij <= dij, etc.

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


function cyclic_triangle_constraints2!(A::SparseMatrixCSC{Float64,Int64},D::Matrix{Int8},X::Matrix{Float64},
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
            ABCproj!(X,W,triplet_corrections,alpha,i,j,k,Wij,Wik,Wjk,n,denom)
            ABCproj!(X,W,triplet_corrections,alpha,i,k,j,Wik,Wij,Wjk,n,denom)
            ABCproj!(X,W,triplet_corrections,alpha,j,k,i,Wjk,Wij,Wik,n,denom)
        end
    end
end

end

function ABCproj!(X::Matrix{Float64},
    W::Matrix{Float64},triplet_corrections::Dict{Int64,Float64},alpha::Float64,
    a::Int64,b::Int64,c::Int64,wab::Float64,wac::Float64,wbc::Float64,n::Int64,denom::Float64)

    # Project onto the constraint convex set defined by xab <= xac + xbc
    # where it is NOT necessarily the case that a < b < c

    ijkKey = (a-1)*n^2+(b-1)*n + c
    corr = pop!(triplet_corrections,ijkKey,0.0)
    Xab = X[b,a]
    Xac = X[max(a,c),min(a,c)]
    Xbc = X[max(b,c),min(b,c)]
    # Correction step

    if corr > 0
        X[b,a] = Xab + corr*wbc*wac/denom
        X[max(a,c),min(a,c)] = Xac - corr*wbc*wab/denom
        X[max(b,c),min(b,c)] = Xbc - corr*wac*wab/denom

        Xab = X[b,a]
        Xac = X[max(a,c),min(a,c)]
        Xbc = X[max(b,c),min(b,c)]
    end

    delta = Xab - Xac - Xbc
    if delta > 1e-8
        X[b,a] = Xab - alpha*delta*wbc*wac/denom
        X[max(a,c),min(a,c)] = Xac + alpha*delta*wab*wbc/denom
        X[max(b,c),min(b,c)] = Xbc + alpha*delta*wab*wac/denom
        triplet_corrections[ijkKey] = delta
    end

end


function cyclic_triangle_constraints!(A::SparseMatrixCSC{Float64,Int64},D::Matrix{Int8},X::Matrix{Float64},Y::Matrix{Float64},
    W::Matrix{Float64},triplet_corrections::Dict{Int64,Float64},P::Matrix{Float64},Q::Matrix{Float64})

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

            ## Check the i,j constraint
            a = i
            b = j
            cor = P[b,a]
            Xab = X[b,a] - cor
            Yab = Y[b,a] - cor
            delta = -Xab - Yab + D[b,a]

            if delta > 0
              X[b,a] = Xab + delta/2
              Y[b,a] = Yab + delta/2
              P[b,a] = delta/2
            else
              P[b,a] = 0.0
            end

            #  X - Y <= -D
            cor = Q[b,a]
            Xab = X[b,a] + cor
            Yab = Y[b,a] - cor
            delta = Xab - Yab - D[b,a]
            if delta > 0
              X[b,a] = Xab - delta/2
              Y[b,a] = Yab + delta/2
              Q[b,a] = delta/2
            else
              Q[b,a] = 0.0
            end

            ## Check the i,j constraint
            a = i
            b = k
            cor = P[b,a]
            Xab = X[b,a] - cor
            Yab = Y[b,a] - cor
            delta = -Xab - Yab + D[b,a]

            if delta > 0
              X[b,a] = Xab + delta/2
              Y[b,a] = Yab + delta/2
              P[b,a] = delta/2
            else
              P[b,a] = 0.0
            end

            #  X - Y <= -D
            cor = Q[b,a]
            Xab = X[b,a] + cor
            Yab = Y[b,a] - cor
            delta = Xab - Yab - D[b,a]
            if delta > 0
              X[b,a] = Xab - delta/2
              Y[b,a] = Yab + delta/2
              Q[b,a] = delta/2
            else
              Q[b,a] = 0.0
            end

            ## Check the i,j constraint
            a = j
            b = k
            cor = P[b,a]
            Xab = X[b,a] - cor
            Yab = Y[b,a] - cor
            delta = -Xab - Yab + D[b,a]

            if delta > 0
              X[b,a] = Xab + delta/2
              Y[b,a] = Yab + delta/2
              P[b,a] = delta/2
            else
              P[b,a] = 0.0
            end

            #  X - Y <= -D
            cor = Q[b,a]
            Xab = X[b,a] + cor
            Yab = Y[b,a] - cor
            delta = Xab - Yab - D[b,a]
            if delta > 0
              X[b,a] = Xab - delta/2
              Y[b,a] = Yab + delta/2
              Q[b,a] = delta/2
            else
              Q[b,a] = 0.0
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
