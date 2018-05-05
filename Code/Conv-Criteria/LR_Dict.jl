#
# LR using new primal-dual checks.
#
# March 21

include("Tricon_Helper.jl")

function ConstructDandW(A::SparseMatrixCSC{Float64,Int64},lam::Float64)
 n = size(A,1)
  Amat = zeros(Float64,n,n)
  W = lam*ones(n,n)
  for i = 1:n-1
      for j = i+1:n
          if A[i,j] > .9
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
        X = -gam*Amat[:,:]
        #W = tril(ones(n,n))
        for i = 1:n
            W[i,i] = 0.0
        end
        # # Matrices for double-loop correction terms (purposely dense)
        P = zeros(n,n)
        SumCorrection1 = 0.0
        SumCorrection2 = 0.0
        InvSumWeights = 0.0
        SlackMax = 0.0
        TotalSlack = 0.0
        for i = 1:n-1
            for j = i+1:n
                InvSumWeights += 1/W[j,i]
            end
        end
        # Correction term vector for triangle constraints
        current_corrections = Vector{Tuple{Int64,Float64}}()
        # next_corrections = Vector{Tuple{Int64,Float64}}()
        # push!(current_corrections,(0,0.0))

        triplet_corrections = Dict{Int64,Float64}()
        # First time through constraints

        cyclic_triangle_constraints!(A,D,X,W,triplet_corrections,1.0)
        #cyclic_triple_loop!(X,W,current_corrections,next_corrections)
        SumCorrection1 = sum_constraint!(X,W,SumCorrection1,InvSumWeights)
        SumCorrection2 = sum_constraint2!(X,W,SumCorrection2,InvSumWeights)
        box_constraints!(X,W,P)


        iter = 0
        lastTri = 1.0
        objective = 0.0
        while iter < maxits
                iter += 1
                vecStart = vec(X)

                # An empty vector to fill future correction terms in
                # current_corrections = next_corrections
                # next_corrections = Vector{Tuple{Int64,Float64}}()
                # push!(current_corrections,(0,0.0))

                tic()
                #cyclic_triple_loop!(X,W,current_corrections,next_corrections)
                cyclic_triangle_constraints!(A,D,X,W,triplet_corrections,1.0)
                SumCorrection1 = sum_constraint!(X,W,SumCorrection1,InvSumWeights)
                SumCorrection2 = sum_constraint2!(X,W,SumCorrection2,InvSumWeights)
                box_constraints!(X,W,P)

                TriTime = toq()
                vecNext = vec(X)
                vecChange = sqrt(mapreduce(z -> abs(z[1]-z[2])^2, +, zip(vecStart,vecNext)))

                # tricheck, lastTri, objective, SlackMax, TotalSlack, SumCorrection =
                # new_report_progress(A,X,P,current_corrections,filename,vecChange,iter,TriTol,lam,lastTri,SlackMax,TotalSlack,TriTime,SumCorrection)

                tricheck, lastTri,objective = report_progress(A,X,W,P,current_corrections,filename,vecChange,iter,TriTol,lam,lastTri,TriTime,SumCorrection1,SumCorrection2,gam)

                if vecChange < Xtol && tricheck && iter > 100
                  break
                end
        end
        Finaltr = FullTriangleCheck(X)
        Finalobj = LR_obj(A,X)
        Finalits = iter
        return X, Finaltr, Finalobj, Finalits
end

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

function LR_obj(A::SparseMatrixCSC{Float64,Int64},X::Matrix{Float64})
  n = size(A,1)
  lr = sum((A[j,i]*X[j,i] for i=1:n-1 for j = i+1:n))
  return lr
end

function report_progress(A::SparseMatrixCSC{Float64,Int64},X::Matrix{Float64},W::Matrix{Float64},P::Matrix{Float64},current_corrections::Vector{Tuple{Int64,Float64}},
    filename::String,xChange::Float64,iter::Int64,TriTol::Float64,lam::Float64,lastTri::Float64,triTime::Float64,SumCorrection1::Float64,SumCorrection2::Float64,gam::Float64)

        n = size(X,1)
        tricheck = TriangleCheck(X,TriTol)
        objective = LR_obj(A,X)
        ch = round(xChange,3)
        ob = round(objective,3)
        time = round(triTime,3)d

        if iter%20 == 1
          # Find the worst constraint violation
          lastTri = FullTriangleCheck(X)

          # Non-metric constraints violations
          sumX = 0.0
          worstnegative = 0.0
          for i = 1:n-1
              for j = i+1:n
                  sumX += X[j,i]
                  if X[j,i] < worstnegative
                      worstnegative = negativity
                  end
              end
          end
          sumVio = abs(n-sumX)

          lastCon = max([lastTri,sumVio,abs(worstnegative)])
        end
        tr = round(lastCon,5)
        nnzDelta = length(current_corrections)

        xWx = xWnorm(W,X)
        BtDual = n*SumCorrection1 - n*SumCorrection2
        PrimalQP = xWx/(2*gam) + LR_obj(A,X)
        DualQP = -BtDual/(gam) - xWx/(2*gam)

        println("$iter: \t Dual = $DualQP, Primal = $PrimalQP, ConSat = $tr, 3Loop = $time, Obj: $ob")

        open(filename, "a") do f
          write(f,"$iter: \t Dual = $DualQP, Primal = $PrimalQP, ConSat = $tr, 3Loop = $time, Obj: $ob\n")
        end

        return tricheck, lastCon, objective
end


function sum_constraint!(X::Matrix{Float64},W::Matrix{Float64},SumCorrection::Float64,InvSumW::Float64)

    # Check if the sum is less than n

    n = size(X,1)
    # Correction step
    sumX = 0.0
    for i = 1:n-1
        for j = i+1:n
            X[j,i] = X[j,i] + SumCorrection/W[j,i]
            sumX += X[j,i]
        end
    end

    if sumX > n
        constant = (sumX - n)/InvSumW
        for i = 1:n-1
            for j = i+1:n
                X[j,i] = X[j,i] - constant/W[j,i]
            end
        end
        return constant
    else
        return 0.0
    end

end


function sum_constraint2!(X::Matrix{Float64},W::Matrix{Float64},SumCorrection::Float64,InvSumW::Float64)

    n = size(X,1)
    # Correction step
    sumX = 0.0
    for i = 1:n-1
        for j = i+1:n
            X[j,i] = X[j,i] - SumCorrection/W[j,i]
            sumX += X[j,i]
        end
    end

    if n > sumX
        constant = (n - sumX)/InvSumW
        for i = 1:n-1
            for j = i+1:n
                X[j,i] = X[j,i] + constant/W[j,i]
            end
        end
        return constant
    else
        return 0.0
    end

end

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

function box_constraints2!(X::Matrix{Float64},W::Matrix{Float64},P::Matrix{Float64})

    n = size(X,1)
    for i = 1:n-1
        for j = i+1:n
            X[j,i] += P[j,i]

            Xij = X[j,i]
            if Xij < 0
                X[j,i] = 0
                P[j,i] = Xij
            else
                P[j,i] = 0.0
            end

        end
    end
end

function cyclic_triangle_constraints!(A::SparseMatrixCSC{Float64,Int64},D::Matrix{Float64},X::Matrix{Float64},
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

# There's a bug in here!!! I have no idea where though.
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
