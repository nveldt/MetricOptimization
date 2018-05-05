# Important note about this code:
# I realized that I think I might be doing something wrong.
# Previously, in the main triple loop (where we check previous corrections,
# and prepare for future projections), in the case where some ordered triplet
# had a nonzero correction but DIDN'T then violate the constraint, nothing changed.
# The actual variabled E[j,i], E[k,i], and E[k,j] only ever changed if, after the
# correction step, they violated their constraint.
#
# This is NOT what is supposed to happen. If there is a nonzero correction,
# that nonzero correction should be applied to really change
# E[j,i], E[k,i], and E[k,j].


# Code for the triangle fixing algorithm of Dhillon et al. specifically
# implemented for LambdaCC
# This version does projections a little differently, as it projects onto a
# weighted 2-norm

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

# FirstTripleLoop!
# The first triple loop to go through on in the Triangle Fixing algorithm
# This is different from other triple loops because it doesn't have any
# adjustment terms.
function FirstTripleLoop!(D::Matrix{Int8},E::Matrix{Float64},W::Matrix{Float64},corrections::Vector{Tuple{Int64,Int64,Int64,Float64}})
  n = size(D,1)
  epsi = 1e-8       # constraint tolerance

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
function TripleLoop2!(D::Matrix{Int8},E::Matrix{Float64},F::Matrix{Float64},P::Matrix{Float64},Q::Matrix{Float64},
  W::Matrix{Float64},now_corrections::Vector{Tuple{Int64,Int64,Int64,Float64}},next_corrections::Vector{Tuple{Int64,Int64,Int64,Float64}})
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
     a = i
     b = j

     cor = P[b,a]
     Eab = E[b,a]- cor
     Fab = F[b,a] - cor
     delta = - Eab - Fab

     if delta > epsi
       E[b,a] = Eab + delta/2
       F[b,a] = Fab + delta/2
       P[b,a] = delta/2
     else
       P[b,a] = 0.0
     end

     #  E - F <= 0

     # corrections
     cor = Q[b,a]
     Eab = E[b,a] - cor
     Fab = F[b,a] + cor
     delta = Eab - Fab
     if delta > epsi
       E[b,a] = Eab - delta/2
       F[b,a] = Fab + delta/2
       Q[b,a] = -delta/2
     else
       Q[b,a] = 0.0
     end

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

    ## Now check double variable constraints

    # i,j
    # a = i
    # b = j
    #
    # cor = P[b,a]
    # Eab = E[b,a]- cor
    # Fab = F[b,a] - cor
    # delta = - Eab - Fab
    #
    # if delta > epsi
    #   E[b,a] = Eab + delta/2
    #   F[b,a] = Fab + delta/2
    #   P[b,a] = delta/2
    # else
    #   P[b,a] = 0.0
    # end
    #
    # #  E - F <= 0
    #
    # # corrections
    # cor = Q[b,a]
    # Eab = E[b,a] - cor
    # Fab = F[b,a] + cor
    # delta = Eab - Fab
    # if delta > epsi
    #   E[b,a] = Eab - delta/2
    #   F[b,a] = Fab + delta/2
    #   Q[b,a] = -delta/2
    # else
    #   Q[b,a] = 0.0
    # end
    #
    # # i,k
    # a = i
    # b = k
    #
    # cor = P[b,a]
    # Eab = E[b,a]- cor
    # Fab = F[b,a] - cor
    # delta = - Eab - Fab
    #
    # if delta > epsi
    #   E[b,a] = Eab + delta/2
    #   F[b,a] = Fab + delta/2
    #   P[b,a] = delta/2
    # else
    #   P[b,a] = 0.0
    # end
    #
    # #  E - F <= 0
    #
    # # corrections
    # cor = Q[b,a]
    # Eab = E[b,a] - cor
    # Fab = F[b,a] + cor
    # delta = Eab - Fab
    # if delta > epsi
    #   E[b,a] = Eab - delta/2
    #   F[b,a] = Fab + delta/2
    #   Q[b,a] = -delta/2
    # else
    #   Q[b,a] = 0.0
    # end
    #
    # # j,k
    # a = j
    # b = k
    #
    # cor = P[b,a]
    # Eab = E[b,a]- cor
    # Fab = F[b,a] - cor
    # delta = - Eab - Fab
    #
    # if delta > epsi
    #   E[b,a] = Eab + delta/2
    #   F[b,a] = Fab + delta/2
    #   P[b,a] = delta/2
    # else
    #   P[b,a] = 0.0
    # end
    #
    # #  E - F <= 0
    #
    # # corrections
    # cor = Q[b,a]
    # Eab = E[b,a] - cor
    # Fab = F[b,a] + cor
    # delta = Eab - Fab
    # if delta > epsi
    #   E[b,a] = Eab - delta/2
    #   F[b,a] = Fab + delta/2
    #   Q[b,a] = -delta/2
    # else
    #   Q[b,a] = 0.0
    # end

end
        a = i
        b = j

        cor = P[b,a]
        Eab = E[b,a]- cor
        Fab = F[b,a] - cor
        delta = - Eab - Fab

        if delta > epsi
          E[b,a] = Eab + delta/2
          F[b,a] = Fab + delta/2
          P[b,a] = delta/2
        else
          P[b,a] = 0.0
        end

        #  E - F <= 0

        # corrections
        cor = Q[b,a]
        Eab = E[b,a] - cor
        Fab = F[b,a] + cor
        delta = Eab - Fab
        if delta > epsi
          E[b,a] = Eab - delta/2
          F[b,a] = Fab + delta/2
          Q[b,a] = -delta/2
        else
          Q[b,a] = 0.0
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
function TriangleFix(A::SparseMatrixCSC{Float64,Int64},D::Matrix{Int8},E::Matrix{Float64},F::Matrix{Float64},W::Matrix{Float64},P::Matrix{Float64},Q::Matrix{Float64},gamma::Float64,Etol::Float64,TriTol::Float64,lam::Float64,filename::String,gam::Float64)
n = size(A,1)
@show filename
  open(filename, "w") do f
          write(f, "Writing output from calling the Triangle Fixing Algorithm\n")
          write(f, "Lambda = $lam, gamma = $gam, Etol = $Etol, TriTol = $TriTol \n")
  end

# Tolerance for inequalty constraints
epsi = 1e-12
tr = 0

# We will fill this with triplets that need corrections in the first round
current_corrections = Vector{Tuple{Int64,Int64,Int64,Float64}}()

# Part 1: Triangle inequality constraint checking
@time FirstTripleLoop!(D,E,W,current_corrections)

# Part 2: -E <= F and E <= F checking, the second set of constraints
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
#Done with the first time through the constraints


# In future rounds, we make corrections with now_corrections, while
# we prepare for new corrections with next_corrections
next_corrections = Vector{Tuple{Int64,Int64,Int64,Float64}}()

iter = 0
maxits = 1e10
while true && iter < maxits
    tempE = E[:,:]

    # Start again with the triangle inequality constraints
    tic()
    TripleLoop2!(D,E,F,P,Q,W,current_corrections,next_corrections)
    lasttime = toq()

    # Get rid of old corrections, prepare space for new
    current_corrections = next_corrections
    TotalViolations = length(current_corrections)
    next_corrections = Vector{Tuple{Int64,Int64,Int64,Float64}}()

    ### Now we check the non-triangle constraints
    # Part 2: -E <= F and E <= F checking

    #DoubleLoop!(P,Q,E,F)


    # Print and update of results for this loop
    iter += 1
    change = norm(vec(E-tempE))
    tricheck = TriangleCheck(E+D,TriTol)
    objective = LPcc_obj(A,E+D,lam)
    ch = round(change,3)
    lt = round(lasttime,1)
    ob = round(objective,3)
    if iter%20 == 1
      G = E+D
      tr = FullTriangleCheck(G)
    end
    println("Iteration $iter, E changed by $ch, 3Loop: $lt, Obj: $ob, Last TriCheck = $tr, Numvio = $TotalViolations")

    open(filename, "a") do f
      write(f,"Iteration $iter, E changed by $ch, 3Loop: $lt, Obj: $ob, Last TriCheck = $tr\n")
    end

    if change < Etol && tricheck
      break
    end

end #end while loop

correctionsLength = length(current_corrections)
next = 1
nextTriplet = current_corrections[next]
tic()
for i = 1:n-2
    for j = i+1:n-2
        for k = j+1:n

            if i == nextTriplet[1] && j == nextTriplet[2] && k == nextTriplet[3]
            cor = current_corrections[next]
              if next < correctionsLength
              next +=1
              nextTriplet = current_corrections[next]
              end
            end

            if i == nextTriplet[1] && k == nextTriplet[2] && j == nextTriplet[3]
            cor = current_corrections[next]
              if next < correctionsLength
              next +=1
              nextTriplet = current_corrections[next]
              end
            end

            if j == nextTriplet[1] && k == nextTriplet[2] && i == nextTriplet[3]
            cor = current_corrections[next]
              if next < correctionsLength
              next +=1
              nextTriplet = current_corrections[next]
              end
            end

        end
    end
end
triloop = toc()

@show triloop

return E+D


end