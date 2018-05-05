include("Tricon_Helper.jl")

# Unless there is a bug, which I don't think there is...this is a huge failure.
# Perhaps because we are not really solving the sub-projections that well, we
# shouldn't hope to solve the overall projections either.
#
# I really don't know why this doesn't work. I would have expected this to work
# at least on tiny graphs like Karate, but perhaps because the convex sets are
# more complicated than just half spaces, convergence is just way worse.

# FirstTripleLoop!
# The first triple loop to go through on in the Triangle Fixing algorithm
# This is different from other triple loops because it doesn't have any
# adjustment terms.
function FirstTripleLoop!(D::Matrix{Int8},E::Matrix{Float64},W::Matrix{Float64},corrections::Vector{Tuple{Int64,Int64,Int64,Float64}})
  n = size(D,1)
  epsi = 0      # constraint tolerance

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
  epsi = 0

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

# Project onto the double constraint
function Project_DoubleConstraints!(D::Matrix{Int8},E::Matrix{Float64},F::Matrix{Float64},tol::Float64)
  println("Constraint 2: Making fij = |xij - dij| \n")

  epsi = 1e-12
  n = size(D,1)
  P = zeros(n,n)
  Q = zeros(n,n)
  change = 1

  while change > tol
    tempE = E[:,:]
    tempF = F[:,:]
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
  change = norm([vec(E-tempE); vec(F-tempF)])
end
  return E,F

end


# This gives the solution for the L2 metric nearness problem
function Project_TriConstraints!(D::Matrix{Int8},E::Matrix{Float64},W::Matrix{Float64},tol::Float64)

  println("Constraint 1: Satisfying triangle inequalities \n")

  n = size(A,1)
  tr = 0
  # We will fill this with triplets that need corrections in the first round
  current_corrections = Vector{Tuple{Int64,Int64,Int64,Float64}}()
  FirstTripleLoop!(D,E,W,current_corrections)

  # In future rounds, we make corrections with now_corrections, while
  # we prepare for new corrections with next_corrections
  next_corrections = Vector{Tuple{Int64,Int64,Int64,Float64}}()

  # Change this eventually
  TriTol = 1e-12

  iter = 0
  maxits = 1e4
  while true && iter < maxits

      tempE = E[:,:]

      # Start again with the triangle inequality constraints
      tic()
      TripleLoop!(D,E,W,current_corrections,next_corrections)
      lasttime = toq()

      # Get rid of old corrections, prepare space for new
      current_corrections = next_corrections
      next_corrections = Vector{Tuple{Int64,Int64,Int64,Float64}}()

      # Print and update of results for this loop
      iter += 1
      change = norm(vec(E-tempE))

      tricheck = TriangleCheck(E+D,TriTol)
      ch = round(change,4)
      lt = round(lasttime,1)
      if iter%20 == 1
        G = E+D
        tr = FullTriangleCheck(G)
      end
      tri = round(tr,5)
    #  println("\t Inner Iteration $iter, E changed by $ch, 3Loop: $lt, Last TriCheck = $tri")

      if change < tol && tricheck
        break
      end

  end #end while loop

  return E

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
function Dykstra_TriangleFix(A::SparseMatrixCSC{Float64,Int64},D::Matrix{Int8},E::Matrix{Float64},F::Matrix{Float64},W::Matrix{Float64},Etol::Float64,TriTol::Float64,lam::Float64,filename::String,gam::Float64)
n = size(A,1)
@show filename
  open(filename, "w") do f
          write(f, "Writing output from calling the Triangle Fixing Algorithm. Two part Dykstra Style\n")
          write(f, "Lambda = $lam, gamma = $gam, Etol = $Etol, TriTol = $TriTol \n")
  end

# Tolerance for inequalty constraints
epsi = 0

# Convergence for inner triple loop
tol1 = 1e-12
tol2 = 1e-12
tr = 1      # tri-check value

# Part 1: Project onto constraint Te <= b
Eold = copy(E)
E = Project_TriConstraints!(D,copy(E),W,tol1)    # Projection
CorrE_tri = copy(E) - copy(Eold)                 # Store correction for next iteration

objective = LPcc_obj(A,E+D,lam)
println("Best l2 norm approximation has objective $objective")
# At this point we should have found the vector is minimum distance from the
# staring vector E and satisfies triangle inequality constraints.
#
# Note: if we were solving the L2 metric nearness problem we'd be done!

# But we are actually interested in solving L1 metric nearness, so we need to
# satisfy the rest of the inequalities as well

# To double check that there aren't mistakes, I copied lots of vectors over and over again.
# This is inefficient, but I wanted to know for sure what was happening with
# different variables and which places in memory they point to. (I know Julia does things
# a little differently from Matlab).

Eold = copy(E)
Fold = copy(F)
output = Project_DoubleConstraints!(D,copy(E),copy(F),tol2)
E = copy(output[1])
F = copy(output[2])
CorrE_doub = copy(E) - copy(Eold)
CorrF_doub = copy(F) - copy(Fold)
objective = LPcc_obj(A,copy(E)+D,lam)
println("Afterwards: objective $objective ")

iter = 0
maxits = 500
while true && iter < maxits

    tempE = E[:,:]
    tempF = F[:,:]

    # Part 1: Project onto constraint Te <= b
    tic()
    Eold = copy(E) - copy(CorrE_tri)
    E = Project_TriConstraints!(D,copy(Eold),W,tol1)
    CorrE_tri = copy(E) - copy(Eold)
    lasttime = toq()
    objective = LPcc_obj(A,copy(E)+D,lam)
    println("Best l2 norm approximation has objective $objective")

    ### Now we check the non-triangle constraints
    # Part 2: -E <= F and E <= F checking
    Eold = copy(E) - copy(CorrE_doub)
    Fold = copy(F) - copy(CorrF_doub)
    output = Project_DoubleConstraints!(D,copy(Eold),copy(Fold),tol2)
    E = copy(output[1])
    F = copy(output[2])
    CorrE_doub = copy(E) - copy(Eold)
    CorrF_doub = copy(F) - copy(Fold)

    objective = LPcc_obj(A,E+D,lam)
    println("After l1 correction: objective $objective")

    # Print and update of results for this loop
    iter += 1
    change = norm([vec(E-tempE); vec(F-tempF)])
    tricheck = TriangleCheck(E+D,TriTol)
    objective = LPcc_obj(A,E+D,lam)
    ch = round(change,5)
    lt = round(lasttime,1)
    ob = round(objective,3)
    if iter%20 == 1
      G = E+D
      tr = FullTriangleCheck(G)
    end

    println("Main Iteration $iter, E changed by $ch, 3Loop: $lt, Obj: $ob, Last TriCheck = $tr")

    open(filename, "a") do f
      write(f,"Main Iteration $iter, E changed by $ch, 3Loop: $lt, Obj: $ob, Last TriCheck = $tr\n")
    end

    if change < Etol && tricheck && iter > 20
      break
    end

end #end while loop

return E+D

end
