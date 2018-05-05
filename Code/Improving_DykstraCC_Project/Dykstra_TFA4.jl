#
# I've used several versions of these algorithms, and now I want to implement
# the cleanest and fastest version of Dykstra's possible to compare all future
# attempts against.
#
#   Coded by Nate Veldt on February 19, 2018

include("Tricon_Helper.jl")

import Base.hash
hash(x::Integer) = UInt64(x)

function ConstructDandW(A::SparseMatrixCSC{Float64,Int64})
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
  return D, W
end

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

        push!(current_corrections,(0,0.0))
        # First time through constraints
        cyclic_triple_loop!(D,E,W,current_corrections,next_corrections)
        cyclic_double_loop!(E,F,P,Q)

        iter = 0
        lastTri = 1.0
        SlackMax = 0.0
        TotalSlack = 0.0
        objective = 0.0
        while iter < maxits
                iter += 1
                vecStart = [vec(E);vec(F)]

                # An empty vector to fill future correction terms in
                current_corrections = next_corrections
                next_corrections = Vector{Tuple{Int64,Float64}}()

                #push!(current_corrections,(0,0,0,0.0))

                tic()
                cyclic_triple_loop!(D,E,W,current_corrections,next_corrections)
                TriTime = toq()

                cyclic_double_loop!(E,F,P,Q)
                vecNext = [vec(E);vec(F)]
                vecChange = sqrt(mapreduce(z -> abs(z[1]-z[2])^2, +, zip(vecStart,vecNext)))
                #vecChange = norm(xStart-xNext)

                tricheck, lastTri, objective, SlackMax, TotalSlack = report_progress(A,E+D,E,F,P,Q,next_corrections,filename,vecChange,iter,TriTol,lam,lastTri,SlackMax,TotalSlack,TriTime)

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


# Check the worst violation in the KKT conditions
function CheckKKT(X::Matrix{Float64},E::Matrix{Float64},F::Matrix{Float64},now_corrections::Vector{Tuple{Int64,Float64}},P::Matrix{Float64},Q::Matrix{Float64})

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
    #
    for i = 1:n-1
        for j = i+1:n
            Eij = E[j,i]
            Fij = F[j,i]
            if Eij - Fij > maxi
                maxi = Eij - Fij
            elseif -Eij - Fij > maxi
                maxi = -Eij - Fij
            end
            slackness = abs((Eij - Fij)*Q[j,i])
            TotalSlack += slackness
            if slackness > SlackMax
                SlackMax = slackness
            end
            slackness = abs((-Eij - Fij)*P[j,i])
            TotalSlack += slackness
            if slackness > SlackMax
                SlackMax = slackness
            end
        end
    end
    # @show SlackMax, TotalSlack
    return round(maxi,3), round(SlackMax,3), round(TotalSlack,3)
  end



function report_progress(A::SparseMatrixCSC{Float64,Int64},X::Matrix{Float64},E::Matrix{Float64},F::Matrix{Float64},P::Matrix{Float64},Q::Matrix{Float64},current_corrections::Vector{Tuple{Int64,Float64}},
    filename::String,xChange::Float64,iter::Int64,TriTol::Float64,lam::Float64,lastTri::Float64,SlackMax::Float64,TotalSlack::Float64,triTime::Float64)

        tricheck = TriangleCheck(X,TriTol)
        objective = LPcc_obj(A,X,lam)
        ch = round(xChange,3)
        ob = round(objective,3)
        time = round(triTime,3)
        tr = round(lastTri,3)
        if iter%10 == 1
         # lastTri = FullTriangleCheck(X)
          lastTri, SlackMax, TotalSlack = CheckKKT(X,E,F,current_corrections,P,Q)
        end
        nnzDelta = length(current_corrections)
        println("Iteration $iter, ||x_k - x_{k+1}|| = $ch, 3Loop = $time, Obj: $ob, ConVio = $lastTri, SlackVio = $SlackMax , T-Slack = $TotalSlack")

        open(filename, "a") do f
          write(f,"Iteration $iter, ||x_k - x_{k+1}|| = $ch, 3Loop = $time, Obj: $ob, Last TriCheck = $tr\n")
        end

        return tricheck, lastTri, objective, SlackMax, TotalSlack
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

function cyclic_triple_loop!(D::Matrix{Float64},E::Matrix{Float64},W::Matrix{Float64},now_corrections::Vector{Tuple{Int64,Float64}},next_corrections::Vector{Tuple{Int64,Float64}})
  n = size(D,1)
  epsi = 1e-8

  correctionsLength = length(now_corrections)
  nowInd = 1;
  # Grab next triplet in the list
  nextKey = now_corrections[nowInd]

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
    end
    ### Done checking triangle j,k,i

end
end
end # end triple for loop

end # end TripleLoop! function
