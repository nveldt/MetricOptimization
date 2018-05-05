# I want to code up the parallel version of Dykstra/Han's method, even though I
# don't actually run it in parallel. This is just to get the different mechanics
# of the algorithm down and compare its convergence, even though only one
# thread will be running it.

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


# Evaluate the LambdaCC LP relaxed objective, given a distance matrix D
function LPcc_obj(A::SparseMatrixCSC{Float64,Int64},D::Matrix{Float64},lam::Float64)
  n = size(A,1)
  # assert(issymmetric(D))
  # assert(issymmetric(A))
  numedges = countnz(A)/2
  lccBound = sum((A[j,i]-lam)*D[j,i] for i=1:n-1 for j = i+1:n) + lam*(n*(n-1)/2 - numedges)
  return lccBound
end

# FirstIteration!
# The first time through the constraints, where we don't need to adjust because of any variables
function FirstIteration!(D::Matrix{Int8},E::Matrix{Float64},F::Matrix{Float64},W::Matrix{Float64},P::Matrix{Float},Q::Matrix{Float},corrections::Vector{Vector{Tuple{Int64,Int64,Int64,Float64}})
  n = size(D,1)
  epsi = 1e-8       # constraint tolerance

  # number of constraints
  numcon = n*(n-1)^2/2

  # Understanding the 'corrections' list:

  # corrections[i] is of tuples of the form (Int, Int, Int, Float)
  # e.g. (1,j,k, mu)
  #
  # Each corresponds to a violation at the triplet i,j,k where i < j < k
  # There are three constraints for triplet (i,j,k):
  # (1)   eij - eik - ejk <= Dik + Djk - Dij
  # (2)  -eij + eik - ejk <= -Dik + Djk + Dij
  # (3)  -eij - eik + ejk <= Dik - Djk + Dij
  #
  # We know what i,j, and k are, and that i < j < k, but we need to know which
  # of these constraints was violated, so that is what we include in the first
  # Int of the tuple--the index of the violated constraint.
  #
  # In this method, we simply iterate through all the constraints in parallel,
  # and record which ones were violated.
  #
  # Later, we will go through and make appropriate adjustments.

# STEP 1: Visit all triangle constraints in parallel and keep track of which
# ones are violated.

# When passing through the triple loop (in parallel), we do not update the value
# of E ever, we only keep track of constraints that are violated
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
    mu = (Eij - Ejk - Eik - b)

    if mu > epsi
      push!(corrections[i],(1,j,k,mu))
      # no further adjustment will be needed for this triplet if this happened
      continue
    end
    # Done checking triangle i,j,k

    # Check triangle i,k,j
    b = -Dik + Djk + Dij
    mu = (-Eij - Ejk + Eik - b)
    if mu > epsi
      push!(corrections[i],(2,j,k,mu))
      continue
    end
    # Done checking triangle i,k,j

    # Check triangle j,k,i
    b = Dik - Djk + Dij
    mu = (-Eij + Ejk - Eik - b)
    if mu > epsi
      push!(corrections[i],(3,j,k,mu))
      continue
    end
    # Done checking triangle j,k,i

end
end
end # end triple for loop


# STEP 2: Adjust Enew based on updates from triplet constraints


# After we have collected all the things that need to be adjusted (in parallel)
# we must go through the list and make appropriate adjustments

# This counts through which index of corrections[i] we are at
Counters = ones(n)

# Initialize a zero matrix
Enew = zeros(n,n)

# each adjustment will be weighted equally among all constraints
lamt = 1/numcon

# Make adjustements and add them to the new matrix G
# This is done sequentially
for a = 1:n
  lengtha = length(corrections[a])

  for v = 1:lengtha
      Trinum = corrections[a][Counters[a]][1]
      i = a
      j = corrections[a][Counters[a]][2]
      k = corrections[a][Counters[a]][3]
      mu = corrections[a][Counters[a]][4]

      # This means there's a triplet (i,j,k), and Trinum tells you which variable
      # Eij, Eik, Ejk is the "odd one out" in the constraint.
      Wik = W[k,i]
      Wjk = W[k,j]
      Wij = W[j,i]
      Eik = E[k,i]
      Ejk = E[k,j]
      Eij = E[j,i]

      denom = Wij*Wjk + Wik*Wij + Wjk*Wik

      signs = ones(3,1)
      signs(Trinum) = -1;

      # Two of them have an added term, and the other has a term subtracted
      Enew[j,i] += lamt*(Eij + signs[1]*mu*(Wik*Wjk/denom))
      Enew[k,i] += lamt*(Eik + signs[2]*mu*(Wij*Wjk/denom))
      Enew[k,j] += lamt*(Ejk + signs[3]*mu*(Wik*Wij/denom))

      # Update how much "lambda-weight" has been added to each entry
      Elam[j,i] += lamt
      Elam[k,i] += lamt
      Elam[k,j] += lamt
  end
end

## Now we handle the constraints of the form Eij - Fij <= 0,
# -Eij - Fij <= 0.
# We do not consider the constraints violated and projections made by
# considering the triangle inequality consraints.
# Here we again go back to what is currenctly stored in E and F originally.
# We store the updates in G, not in E. And we scale the update by the corresponding
# lambda weight.

Flam = zeros(n,n)
Fnew = zeros(n,n)

for i = 1:n-1
  for j = i+1:n
    Eij = E[j,i]
    Fij = F[j,i]
    delta = -Eij - Fij
    if delta > epsi
      Enew[j,i] += lamt*(Eij + delta/2)
      Fnew[j,i] += lamt*(Fij + delta/2)
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
      Enew[j,i] += lamt*(Eij - delta/2)
      Fnew[j,i] += lamt*(Fij + delta/2)
      Q[j,i] = -delta/2
      # else the correction next time will just be P[i,j] == 0
    end
  end
end

# We have only incremented entry (i,j) of Enew and Fnew
# due to the violated constraints that they were a part of. We now need to
# update Eij with the remaining "unused" lambda-weight times the previous
# value of entry i,j
for i = 1:n-1
  for j = i+1:n
    Enew[j,i] += (1-Elam[j,i])*E[j,i]
    Fnew[j,i] += (1-Flam[j,i])*F[j,i]
  end
end

# Done with the first time through the constraints


end # end FirstTripleLoop! function


# MainIteration!
# The first time through the constraints, where we don't need to adjust because of any variables
function MainIteration!(D::Matrix{Int8},E::Matrix{Float64},F::Matrix{Float64},W::Matrix{Float64},P::Matrix{Float},Q::Matrix{Float},L_curr::Vector{Vector{Tuple{Int64,Int64,Int64,Float64}},L_next::Vector{Vector{Tuple{Int64,Int64,Int64,Float64}})
  n = size(D,1)
  epsi = 1e-8       # constraint tolerance

  # number of constraints
  numcon = n*(n-1)^2/2

  # Understanding the 'L_next' list:

  # L_next[i] is of tuples of the form (Int, Int, Int, Float)
  # e.g. (1,j,k, mu)
  #
  # Each corresponds to a violation at the triplet i,j,k where i < j < k
  # There are three constraints for triplet (i,j,k):
  # (1)   eij - eik - ejk <= Dik + Djk - Dij
  # (2)  -eij + eik - ejk <= -Dik + Djk + Dij
  # (3)  -eij - eik + ejk <= Dik - Djk + Dij
  #
  # We know what i,j, and k are, and that i < j < k, but we need to know which
  # of these constraints was violated, so that is what we include in the first
  # Int of the tuple--the index of the violated constraint.
  #
  # In this method, we simply iterate through all the constraints in parallel,
  # and record which ones were violated.
  #
  # Later, we will go through and make appropriate adjustments.

# STEP 1: Visit all triangle constraints in parallel and keep track of which
# ones are violated.

# We have n different lists L_curr[i] for i = 1,2, ...,n
# Each requires a different counter to we can keep track of which triplet
# we are currently visiting. There are n "nextTriplets" at a time, which we also
# must keep track of


# Initialize the nextTriplets vector, and vector of lengths of lists
ListLength = zeros(n);

# This probably won't work...need to come back later to deal with nextTriplet
nextTriplet = Vector{Tuple{Int64,Int64,Int64,Float64}}()
for i = 1:n
  #push!(nextTriplet[i],L_curr[i][1])
  nextTriplet[i] = L_curr[i][1]
  ListLength[i] = length(L_curr[i])
end

Counters = ones(n)

# When passing through the triple loop (in parallel), we do not update the value
# of E ever, we only keep track of constraints that are violated
#
# Before checking for constraint violation, we adjust the variabled Eij, Eik,
# and Ejk if there is a nonzero correction/adjustment value in L_curr
@inbounds for i = 1:n-2
   for j = i+1:n-1
   Eij = E[j,i]
   Dij = D[j,i]
   Wij  = W[j,i]
    for k = j+1:n
    Dik = D[k,i]
    Djk = D[k,j]
    Eik = E[k,i]
    Ejk = E[k,j]
    Wik = W[k,i]
    Wjk = W[k,j]

    # Consider defining a new function here that you can just pass things off
    # to three times so that you don't have to repeat a bunch of code
    #
    # Need to understand how different threads will treat the fact that we're
    # re-naming Eik, Dik, etc. every time. Will these be local to the thread
    # that defines them?

    # Check triplet i,j,k
    if j == nextTriplet[i][2] && k == nextTriplet[i][3] && 1 == nextTriplet[i][1]
      cor = nextTriplet[i][4]

      # We need to scale the correction since we're projecting into the minimum
      # weighted 2-norm direction
      denom = Wij*Wjk + Wik*Wij + Wjk*Wik

      Eij = Eij + cor*(Wik*Wjk/denom)
      Eik = Eik - cor*(Wij*Wjk/denom)
      Ejk = Ejk - cor*(Wik*Wij/denom)

      # Move along in the list of triplets with corrections
      if Counters[i] < ListLength[i]
        Counters[i] +=1
        nextTriplet[i] = L_curr[i][Counters[i]]
      end
    end

    b = Dik + Djk - Dij
    mu = (Eij - Ejk - Eik - b)

    if mu > epsi
      push!(L_next[i],(1,j,k,mu))
      # no further adjustment will be needed for this triplet if this happened
      continue
    end
    # Done checking triangle i,j,k

    # Check triplet i,k,j
    if j == nextTriplet[i][2] && k == nextTriplet[i][3] && 1 == nextTriplet[i][2]
      cor = nextTriplet[i][4]

      # We need to scale the correction since we're projecting into the minimum
      # weighted 2-norm direction
      denom = Wij*Wjk + Wik*Wij + Wjk*Wik

      Eij = Eij - cor*(Wik*Wjk/denom)
      Eik = Eik + cor*(Wij*Wjk/denom)
      Ejk = Ejk - cor*(Wik*Wij/denom)

      # Move along in the list of triplets with corrections
      if Counters[i] < ListLength[i]
        Counters[i] +=1
        nextTriplet[i] = L_curr[i][Counters[i]]
      end
    end

    b = -Dik + Djk + Dij
    mu = (-Eij - Ejk + Eik - b)

    if mu > epsi
      push!(L_next[i],(2,j,k,mu))
      # no further adjustment will be needed for this triplet if this happened
      continue
    end
    # Done checking triangle i,k,j

    # Check triplet j,k,i
    if j == nextTriplet[i][2] && k == nextTriplet[i][3] && 1 == nextTriplet[i][3]
      cor = nextTriplet[i][4]

      # We need to scale the correction since we're projecting into the minimum
      # weighted 2-norm direction
      denom = Wij*Wjk + Wik*Wij + Wjk*Wik

      Eij = Eij - cor*(Wik*Wjk/denom)
      Eik = Eik - cor*(Wij*Wjk/denom)
      Ejk = Ejk + cor*(Wik*Wij/denom)

      # Move along in the list of triplets with corrections
      if Counters[i] < ListLength[i]
        Counters[i] +=1
        nextTriplet[i] = L_curr[i][Counters[i]]
      end
    end

    b = Dik - Djk + Dij
    mu = (-Eij + Ejk - Eik - b)

    if mu > epsi
      push!(L_next[i],(3,j,k,mu))
      # no further adjustment will be needed for this triplet if this happened
      continue
    end
    # Done checking triangle j,k,i

end
end
end # end triple for loop


# STEP 2: Adjust Enew based on updates from triplet constraints


# After we have collected all the things that need to be adjusted (in parallel)
# we must go through the list and make appropriate adjustments

# This counts through which index of corrections[i] we are at
Counters = ones(n)

# Initialize a zero matrix
Enew = zeros(n,n)

# each adjustment will be weighted equally among all constraints
lamt = 1/numcon

# Make adjustements and add them to the new matrix G
# This is done sequentially
for a = 1:n
  lengtha = length(corrections[a])

  for v = 1:lengtha
      Trinum = corrections[a][Counters[a]][1]
      i = a
      j = corrections[a][Counters[a]][2]
      k = corrections[a][Counters[a]][3]
      mu = corrections[a][Counters[a]][4]

      # This means there's a triplet (i,j,k), and Trinum tells you which variable
      # Eij, Eik, Ejk is the "odd one out" in the constraint.
      Wik = W[k,i]
      Wjk = W[k,j]
      Wij = W[j,i]
      Eik = E[k,i]
      Ejk = E[k,j]
      Eij = E[j,i]

      denom = Wij*Wjk + Wik*Wij + Wjk*Wik

      signs = ones(3,1)
      signs(Trinum) = -1;

      # Two of them have an added term, and the other has a term subtracted
      Enew[j,i] += lamt*(Eij + signs[1]*mu*(Wik*Wjk/denom))
      Enew[k,i] += lamt*(Eik + signs[2]*mu*(Wij*Wjk/denom))
      Enew[k,j] += lamt*(Ejk + signs[3]*mu*(Wik*Wij/denom))

      # Update how much "lambda-weight" has been added to each entry
      Elam[j,i] += lamt
      Elam[k,i] += lamt
      Elam[k,j] += lamt
  end
end

## Now we handle the constraints of the form Eij - Fij <= 0,
# -Eij - Fij <= 0.
# We do not consider the constraints violated and projections made by
# considering the triangle inequality consraints.
# Here we again go back to what is currenctly stored in E and F originally.
# We store the updates in G, not in E. And we scale the update by the corresponding
# lambda weight.

Flam = zeros(n,n)
Fnew = zeros(n,n)

for i = 1:n-1
  for j = i+1:n
    Eij = E[j,i]
    Fij = F[j,i]
    delta = -Eij - Fij
    if delta > epsi
      Enew[j,i] += lamt*(Eij + delta/2)
      Fnew[j,i] += lamt*(Fij + delta/2)
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
      Enew[j,i] += lamt*(Eij - delta/2)
      Fnew[j,i] += lamt*(Fij + delta/2)
      Q[j,i] = -delta/2
      # else the correction next time will just be P[i,j] == 0
    end
  end
end

# We have only incremented entry (i,j) of Enew and Fnew
# due to the violated constraints that they were a part of. We now need to
# update Eij with the remaining "unused" lambda-weight times the previous
# value of entry i,j
for i = 1:n-1
  for j = i+1:n
    Enew[j,i] += (1-Elam[j,i])*E[j,i]
    Fnew[j,i] += (1-Flam[j,i])*F[j,i]
  end
end

# Done with the first time through the constraints


end # end MainIteration

# TripleLoop!
# This is the interior triple loop
#
# now_corrections is a vector storing which triplets (i,j,k) were changed
# in the last iteration, so we have a correction variable we need to now apply
function TripleLoop!(D::Matrix{Int8},E::Matrix{Float64},W::Matrix{Float64},now_corrections::Vector{Tuple{Int64,Int64,Int64,Float64}},next_corrections::Vector{Tuple{Int64,Int64,Int64,Float64}})
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

      Eij = Eij + cor*(Wik*Wjk/denom)
      Eik = Eik - cor*(Wij*Wjk/denom)
      Ejk = Ejk - cor*(Wik*Wij/denom)

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
      Eij = E[j,i] - cor*(Wik*Wjk/denom)
      Eik = E[k,i] + cor*(Wij*Wjk/denom)
      Ejk = E[k,j] - cor*(Wik*Wij/denom)

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

      Eij = E[j,i] - cor*(Wik*Wjk/denom)
      Eik = E[k,i] - cor*(Wij*Wjk/denom)
      Ejk = E[k,j] + cor*(Wik*Wij/denom)

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
end


# TriangleFix
#
# This solves the LambdaCC LP relaxation using the Triangle Fixing algorithm
# of Dhillon et al. for the L1 metric nearness problem.
#
# D = 0,1 matrix indicating node distances--not satisfying triangle inequalities
# W = weights for each edge in the graph (based on what lambda is)
# gamma = parameter for turning the LP into a 2-norm problem (if this is not
#         set to a high enough value, this problem isn't actually equivalent to
#         the LP relaxation of correlation clustering.)
#
# The code here is basically just a wrapper for the FirstIteration and
# MainIteration functions
function TriangleFix(A::SparseMatrixCSC{Float64,Int64},D::Matrix{Int8},E::Matrix{Float64},F::Matrix{Float64},W::Matrix{Float64},P::Matrix{Float64},Q::Matrix{Float64},gamma::Float64,Etol::Float64,TriTol::Float64,lam::Float64,filename::String,gam::Float64)
  n = size(A,1)
  @show filename
  open(filename, "w") do f
          write(f, "Writing output from calling the Triangle Fixing Algorithm\n")
          write(f, "Lambda = $lam, gamma = $gam, Etol = $Etol, TriTol = $TriTol \n")
  end

  # Tolerance for inequalty constraints
  epsi = 1e-8

  # Define a vector of vectors of tuples, so that for each i we have a list
  # of varaible corrections having to do with triplets (i,j,k) where i < j, i<k
  L_current = Vector{Vector{Tuple{Int64,Int64,Int64,Float64}}}()
  for i = 1:n
      push!(L_current,Vector{Tuple{Int64,Int64,Int64,Float64}}())
  end

  # number of constraints
  numcon = n*(n-1)^2/2

  # Go through the first iteration, and update E and F
  # L_current will store triplet constraints that have been adjusted
  # and P and Q will store double constraints that have been adjusted
  @time FirstIteration!(D,E,F,W,P,Q,L_current)

  # In future rounds, L_current stores corrections that need to be made before
  # making projections, while L_next keeps track of new constraint violations
  # we need to keep track of
  L_next = Vector{Vector{Tuple{Int64,Int64,Int64,Float64}}}()
  for i = 1:n
      push!(L_next,Vector{Tuple{Int64,Int64,Int64,Float64}}())
  end

  iter = 0
  maxits = 300
  while true && iter < maxits
      tempE = E[:,:]

      # Start again with the triangle inequality constraints
      tic()

      MainIteration!(D,E,F,W,P,Q,L_next,L_current)

      lasttime = toq()

      # Get rid of old corrections, prepare space for new
      L_current = L_next
      L_next = Vector{Vector{Tuple{Int64,Int64,Int64,Float64}}}()
      for i = 1:n
          push!(L_next,Vector{Tuple{Int64,Int64,Int64,Float64}}())
      end

      # Print an update of results for this loop
      iter += 1
      change = norm(vec(E-tempE))
      tricheck = TriangleCheck(E+D,TriTol)
      objective = LPcc_obj(A,E+D,lam)
      ch = round(change,3)
      lt = round(lasttime,1)
      ob = round(objective,3)
      if iter%20 == 1
        G = E+D
      @time tr = FullTriangleCheck(G)
      end
      println("Iteration $iter, E changed by $ch, 3Loop: $lt, Obj: $ob")

      open(filename, "a") do f
        write(f,"Iteration $iter, E changed by $ch, 3Loop: $lt, Obj: $ob, Last TriCheck = $tr\n")
      end

      if change < Etol && tricheck
        break
      end

  end #end while loop

  return E+D
end
