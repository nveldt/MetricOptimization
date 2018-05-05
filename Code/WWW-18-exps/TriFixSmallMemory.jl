# Code for the triangle fixing algorithm of Dhillon et al. specifically
# implemented for LambdaCC

# This keeps track of corrections when we are iterating through the constraints
# that need to be fixed along the way

include("TriFixHelper.jl")

#
# D = 0,1 matrix indicating node distances--not satisfying triangle inequalities
# W = weights for each edge in the graph
# gamma = parameter for turning the LP into a 2-norm problem
#
function TriangleFix(A::SparseMatrixCSC{Float64,Int64},D::Matrix{Int8},E::Matrix{Float64},F::Matrix{Float64},W::Matrix{Float64},P::Matrix{Float64},Q::Matrix{Float64},gamma::Float64,Etol::Float64,TriTol::Float64,lam::Float64,filename::String,gam::Float64)
n = size(A,1)
  open(filename, "w") do f
          write(f, "Writing output from calling the Triangle Fixing Algorithm\n")
          write(f, "Lambda = $lam, gamma = $gam, Etol = $Etol, TriTol = $TriTol \n")
  end

epsi = 1e-8

# We will fill this with triplets that need corrections in the first round
current_corrections = Vector{Tuple{Int64,Int64,Int64,Float64}}()

# Part 1: Triangle inequality constraint checking
@time FirstTripleLoop!(D,E,current_corrections)

# Part 2: -E <= F and E <= F checking
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
# Done with the first time through the constraints


# In future rounds, we make corrections with now_corrections, while
# we prepare for new corrections with next_corrections
next_corrections = Vector{Tuple{Int64,Int64,Int64,Float64}}()

iter = 0
while true

    tempE = E[:,:]

    # Start again with the triangle inequality constraints
    tic()
    TripleLoop!(D,E,current_corrections,next_corrections)
    lasttime = toq()

    # Get rid of old corrections, prepare space for new
    current_corrections = next_corrections
    next_corrections = Vector{Tuple{Int64,Int64,Int64,Float64}}()


    ### Now we check the non-triangle constraints
    # Part 2: -E <= F and E <= F checking
    DoubleLoop!(P,Q,E,F)


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
