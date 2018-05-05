# Returns whether or not all triangle inequality constraints are satisfied
# to within the desired tolerance
function TriangleCheck(D::Matrix{Float64},tol::Float64)
  # Checks whether or not the triangle inequality constraints for matrix D
  # are satisfied to within the given tolerance
  n = size(D,1)
  answer = true;
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


# Evaluate the LambdaCC relaxed objective

function LPcc_obj(A::SparseMatrixCSC{Float64,Int64},D::Matrix{Float64},lam::Float64)
n = size(A,1)
# assert(issymmetric(D))
# assert(issymmetric(A))
numedges = countnz(A)/2
lccBound = sum((A[j,i]-lam)*D[j,i] for i=1:n-1 for j = i+1:n) + lam*(n*(n-1)/2 - numedges)
return lccBound
end


# The first triple loop to go through on in the Triangle Fixing algorithm
function FirstTripleLoop!(D::Matrix{Int8},E::Matrix{Float64},corrections::Vector{Tuple{Int64,Int64,Int64,Float64}})
n = size(D,1)
epsi = 1e-8

@inbounds for i = 1:n-2
   for j = i+1:n-1
   Eij = E[j,i]
   Dij = D[j,i]
    for k = j+1:n

    Dik = D[k,i]
    Djk = D[k,j]
    Eik = E[k,i]
    Ejk = E[k,j]

    # Check triangle i,j,k
    b = Dik + Djk - Dij
    mu = (Eij - Ejk - Eik - b)/3

    if mu > epsi
      E[j,i] = Eij - mu
      E[k,i] = Eik + mu
      E[k,j] = Ejk + mu

      # Next time we see this triple we have to correct
      push!(corrections,(i,j,k,mu))
      # no further adjustment will be needed for this triplet if this happened
      continue
    end
    # Done checking triangle i,j,k

    # Note: we don't need to reset Eij = E[j,i] etc now. If mu > 0 above, then
    # it is easy to show that no adjustment will be needed for the next two
    # orderings of the triplet i,j,k. If mu <= 0, then we didn't update E, so
    # there is no need to reset Eij = E[j,i] etc.


    # Check triangle i,k,j
    b = -Dik + Djk + Dij
    mu = (-Eij - Ejk + Eik - b)/3
    if mu > epsi
      E[j,i] = Eij + mu
      E[k,i] = Eik - mu
      E[k,j] = Ejk + mu

      # Next time we see this triple we have to correct
      push!(corrections,(i,k,j,mu))
      continue
    end
    # Done checking triangle i,k,j

    # Check triangle j,k,i
    b = Dik - Djk + Dij
    mu = (-Eij + Ejk - Eik - b)/3
    if mu > epsi
      E[j,i] = Eij + mu
      E[k,i] = Eik + mu
      E[k,j] = Ejk - mu

      # Next time we see this triple we have to correct
      push!(corrections,(j,k,i,mu))
      continue
    end
    # Done checking triangle j,k,i

end
end
end # end triple for loop

end # end FirstTripleLoop! function


function TripleLoop!(D::Matrix{Int8},E::Matrix{Float64},now_corrections::Vector{Tuple{Int64,Int64,Int64,Float64}},next_corrections::Vector{Tuple{Int64,Int64,Int64,Float64}})
n = size(D,1)
epsi = 1e-8

correctionsLength = length(now_corrections)
nowInd = 1;
# Grab next triplet in the list

if correctionsLength > 0
nextTriplet = now_corrections[nowInd]
else
  nextTriplet = [0, 0 , 0]
end
#nextTriplet = shift!(now_corrections)

@inbounds for i = 1:n-2
  for  j = i+1:n-1
   Eij = E[j,i]
   Dij = D[j,i]
  for  k = j+1:n

    Dik = D[k,i]
    Djk = D[k,j]
    Eik = E[k,i]
    Ejk = E[k,j]

    ### Check triangle i,j,k

    # First see if this is the next triangle with a nonzero correction variable
    if i == nextTriplet[1] && j == nextTriplet[2] && k == nextTriplet[3]
      cor = nextTriplet[4]
      Eij = Eij + cor
      Eik = Eik - cor
      Ejk = Ejk - cor

      # Move along in the list of triplets with corrections
      if nowInd < correctionsLength
      nowInd +=1
      nextTriplet = now_corrections[nowInd]
      end
    end

    b = Dik + Djk - Dij
    mu = (Eij - Ejk - Eik - b)/3

    if mu > epsi
      E[j,i] = Eij - mu
      E[k,i] = Eik + mu
      E[k,j] = Ejk + mu

      # Next time we see this triple we have to correct
      push!(next_corrections,(i,j,k,mu))
    end
    ### Done checking triangle i,j,k


    ### Check triangle i,k,j
    if i == nextTriplet[1] && k == nextTriplet[2] && j == nextTriplet[3]
      cor = nextTriplet[4]
      Eij = E[j,i] - cor
      Eik = E[k,i] + cor
      Ejk = E[k,j] - cor

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
    mu = (-Eij - Ejk + Eik - b)/3
    if mu > epsi
      E[j,i] = Eij + mu
      E[k,i] = Eik - mu
      E[k,j] = Ejk + mu

      # Next time we see this triple we have to correct
      push!(next_corrections,(i,k,j,mu))
    end
    ### Done checking triangle i,k,j

    ### Check triangle j,k,i
    if j == nextTriplet[1] && k == nextTriplet[2] && i == nextTriplet[3]
      cor = nextTriplet[4]
      Eij = E[j,i] - cor
      Eik = E[k,i] - cor
      Ejk = E[k,j] + cor

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
    mu = (-Eij + Ejk - Eik - b)/3
    if mu > epsi
      E[j,i] = Eij + mu
      E[k,i] = Eik + mu
      E[k,j] = Ejk - mu

      # Next time we see this triple we have to correct
      push!(next_corrections,(j,k,i,mu))
    end
    ### Done checking triangle j,k,i

end
end
end # end triple for loop

end # end TripleLoop! function

# Double loop algorithm
function DoubleLoop!(P::Matrix{Float64},Q::Matrix{Float64},E::Matrix{Float64},F::Matrix{Float64})
epsi = 1e-8
n = size(P,1)
@inbounds for i = 1:n-1
  for j = i+1:n

    #  -E - F <= 0

    # corrections
    cor = P[j,i]
    Eij = E[j,i] - cor
    Fij = F[j,i] - cor
    delta = - Eij - Fij

    if delta > epsi
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
    if delta > epsi
      E[j,i] = Eij - delta/2
      F[j,i] = Fij + delta/2
      Q[j,i] = -delta/2
    else
      Q[j,i] = 0.0
    end

  end
end

end
