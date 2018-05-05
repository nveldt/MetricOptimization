function ConstructDandW(A::SparseMatrixCSC{Float64,Int64})
  D = zeros(Int8,n,n)
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

# A Helper Library so that I can keep my other .jl files a little cleaner
# TriangleCheck
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

# FullTriangleCheck
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


# Evaluate the General Correlation Clustering objective, where W is an affinity matrix (has positive and negative edges)
function GeneralCCLP_obj(W::Matrix{Float64},D::Matrix{Float64})
  n = size(W,1)
  LPbound = 0.0
    for i = 1:n-1
        for j = i+1:n
            if W[i,j] > 0
                LPbound += W[i,j]*D[i,j]
            else
                LPbound += -W[i,j]*(1-D[i,j])
            end
        end
    end

  return LPbound
end


# Evaluate the LambdaCC LP relaxed objective, given a distance matrix D
function LPcc_obj(A::SparseMatrixCSC{Float64,Int64},D::Matrix{Float64},lam::Float64)
  n = size(A,1)
  # assert(issymmetric(D))
  # assert(issymmetric(A))
  numedges = countnz(A)/2
  lccBound = sum((A[j,i]-lam)*D[j,i] for i=1:n-1 for j = i+1:n) + lam*(n*(n-1)/2 - numedges)
  return lccBound
end

function LR_obj(A::SparseMatrixCSC{Float64,Int64},X::Matrix{Float64})
  n = size(A,1)
  lr = sum((A[j,i]*X[j,i] for i=1:n-1 for j = i+1:n))
  return lr
end

# This function visits each triangle in turn and projects onto it.
# No correction terms; just project and move on.
function TripleLoop_NoMemory!(D::Matrix{Int8},E::Matrix{Float64},W::Matrix{Float64})
n = size(D,1)
epsi = 0      # constraint tolerance

@inbounds for i = 1:n-2
   for j = i+1:n-1


    for k = j+1:n
    Eij = E[j,i]
    Dij = D[j,i]
    Wij = W[j,i]
    Dik = D[k,i]
    Djk = D[k,j]
    Eik = E[k,i]
    Ejk = E[k,j]
    Wik = W[k,i]
    Wjk = W[k,j]


    # Check triangle i,j,k
    b = Dik + Djk - Dij

    # changed, by removing the 3
    mu = (Eij - Ejk - Eik - b)

    if mu > epsi

      denom = Wij*Wjk + Wik*Wij + Wjk*Wik
      E[j,i] = Eij - mu*(Wik*Wjk/denom)
      E[k,i] = Eik + mu*(Wij*Wjk/denom)
      E[k,j] = Ejk + mu*(Wik*Wij/denom)

    end

    # Done checking triangle i,j,k

    # Note: we don't need to reset Eij = E[j,i] etc now. If mu > 0 above, then
    # it is easy to show that no adjustment will be needed for the next two
    # orderings of the triplet i,j,k. If mu <= 0, then we didn't update E, so
    # there is no need to reset Eij = E[j,i] etc.

    # Check triangle i,k,j
    b = -Dik + Djk + Dij
    mu = (-Eij - Ejk + Eik - b)
    if mu > epsi

      denom = Wij*Wjk + Wik*Wij + Wjk*Wik

      E[j,i] = Eij + mu*(Wik*Wjk/denom)
      E[k,i] = Eik - mu*(Wij*Wjk/denom)
      E[k,j] = Ejk + mu*(Wik*Wij/denom)

    end
    # Done checking triangle i,k,j

    # Check triangle j,k,i
    b = Dik - Djk + Dij
    mu = (-Eij + Ejk - Eik - b)
    if mu > epsi

      denom = Wij*Wjk + Wik*Wij + Wjk*Wik

      E[j,i] = Eij + mu*(Wik*Wjk/denom)
      E[k,i] = Eik + mu*(Wij*Wjk/denom)
      E[k,j] = Ejk - mu*(Wik*Wij/denom)
    end
    # Done checking triangle j,k,i

#    F = -gam*ones(n,n)*sigmaK + (1-sigmaK)*F

end
end
end # end triple for loop

end # end TripleLoop_NoMemory


# Double loop algorithm, no keeping track of previous or future "corrections"
function DoubleLoop_NoMemory!(E::Matrix{Float64},F::Matrix{Float64})
  epsi = 0
  n = size(E,1)
    @inbounds for i = 1:n-1
  for j = i+1:n

    #  -E - F <= 0
    Eij = E[j,i]
    Fij = F[j,i]
    delta = - Eij - Fij

    if delta > epsi
      E[j,i] = Eij + delta/2
      F[j,i] = Fij + delta/2
    end

    #  E - F <= 0
    Eij = E[j,i]
    Fij = F[j,i]
    delta = Eij - Fij
    if delta > epsi
      E[j,i] = Eij - delta/2
      F[j,i] = Fij + delta/2
    end

  end
end

end


# FirstTripleLoop!
# The first triple loop to go through on in the Triangle Fixing algorithm
# This is different from other triple loops because it doesn't have any
# adjustment terms.
function FirstTripleLoop!(D::Matrix{Int8},E::Matrix{Float64},W::Matrix{Float64},corrections::Vector{Tuple{Int64,Int64,Int64,Float64}})
  n = size(D,1)
  epsi = 1e-12       # constraint tolerance

@inbounds for i = 1:n-2
   for j = i+1:n-1
   Eij = E[j,i]
   Dij = D[j,i]
   Wij = W[j,i]
    for k = j+1:n

    Dik = D[k,i]
    Djk = D[k,j]
    Eik = E[k,i]
    Ejk = E[k,j]
    Wik = W[k,i]
    Wjk = W[k,j]

    # Check triangle i,j,k
    b = Dik + Djk - Dij

    # changed, by removing the 3
    mu = (Eij - Ejk - Eik - b)

    if mu > epsi

      denom = Wij*Wjk + Wik*Wij + Wjk*Wik
    #  @show Wik*Wjk/denom
      E[j,i] = Eij - mu*(Wik*Wjk/denom)
      E[k,i] = Eik + mu*(Wij*Wjk/denom)
      E[k,j] = Ejk + mu*(Wik*Wij/denom)

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
    mu = (-Eij - Ejk + Eik - b)
    if mu > epsi

      denom = Wij*Wjk + Wik*Wij + Wjk*Wik

      E[j,i] = Eij + mu*(Wik*Wjk/denom)
      E[k,i] = Eik - mu*(Wij*Wjk/denom)
      E[k,j] = Ejk + mu*(Wik*Wij/denom)

      # Next time we see this triple we have to correct
      push!(corrections,(i,k,j,mu))
      continue
    end
    # Done checking triangle i,k,j

    # Check triangle j,k,i
    b = Dik - Djk + Dij
    mu = (-Eij + Ejk - Eik - b)
    if mu > epsi

      denom = Wij*Wjk + Wik*Wij + Wjk*Wik

      E[j,i] = Eij + mu*(Wik*Wjk/denom)
      E[k,i] = Eik + mu*(Wij*Wjk/denom)
      E[k,j] = Ejk - mu*(Wik*Wij/denom)

      # Next time we see this triple we have to correct
      push!(corrections,(j,k,i,mu))
      continue
    end
    # Done checking triangle j,k,i

end
end
end # end triple for loop

end # end FirstTripleLoop! function


# TripleLoop!
# This is the interior triple loop
#
# now_corrections is a vector storing which triplets (i,j,k) were changed
# in the last iteration, so we have a correction variable we need to now apply
function TripleLoop!(D::Matrix{Int8},E::Matrix{Float64},W::Matrix{Float64},now_corrections::Vector{Tuple{Int64,Int64,Int64,Float64}},next_corrections::Vector{Tuple{Int64,Int64,Int64,Float64}})
  n = size(D,1)
  epsi = 1e-12

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

      # We need to scale the correction since we're projecting into the minimum
      # weighted 2 norm direction
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




# Double loop algorithm
function DoubleLoop!(P::Matrix{Float64},Q::Matrix{Float64},E::Matrix{Float64},F::Matrix{Float64})
  epsi = 1e-12
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

end # End DoubleLoop function

function FirstDoubleLoop!(P::Matrix{Float64},Q::Matrix{Float64},E::Matrix{Float64},F::Matrix{Float64})
  epsi = 0
  n = size(E,1)
  for i = 1:n-1
    for j = i+1:n

      Eij = E[j,i]
      Fij = F[j,i]
      delta = -Eij - Fij
      if delta > epsi
        E[j,i] = Eij + delta/2
        F[j,i] = Fij + delta/2
        P[j,i] = delta/2
        # else the correction next time will just be P[i,j] == 0
      end
  end
  end

  for i = 1:n-1
    for j = i+1:n
      Eij = E[j,i]
      Fij = F[j,i]
      delta = Eij - Fij
      if delta > epsi
        E[j,i] = Eij - delta/2
        F[j,i] = Fij + delta/2
        Q[j,i] = -delta/2
        # else the correction next time will just be P[i,j] == 0
      end
    end
  end
end

function Bauschke_Update!(E::Matrix{Float64},F::Matrix{Float64},sigmaK::Float64,n::Int64,gam::Float64)
  for i = 1:n-1
    for j = i+1:n
      E[j,i] = (1-sigmaK)*E[j,i]
      F[j,i] = -gam*sigmaK + (1-sigmaK)*F[j,i]
    end
  end
end
