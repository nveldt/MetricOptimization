# Go through the triple loop the first time, and identify
# triangle inequality violations that must be fixed
function FirstTripleLoop!(E::Matrix{Float64},D::Matrix{Int8},corrections::Vector{Vector{Tuple{Int64,Int64,Int64,Float64}}})
    n = size(E,1)
     Threads.@threads for i = 1:n-2
        JobI!(E,D,i,corrections[i])
    end
end

function JobI!(E::Matrix{Float64},D::Matrix{Int8},i::Int64,local_corr::Vector{Tuple{Int64,Int64,Int64,Float64}})
    # This is a task performed by a single thread. It is the inner two loops
    # of the triple for loop
    n = size(E,1); epsi = 1e-8
    for j = i+1:n-1
        Xij = E[j,i] + D[j,i]
         for k = j+1:n
         Xik = E[k,i] + D[k,i]; Xjk = E[k,j] + D[k,j]

         mu = Xij - Xik - Xjk
         if mu > epsi
           push!(local_corr,(1,j,k,mu))
           continue
         end

         mu = -Xij + Xik - Xjk
         if mu > epsi
           push!(local_corr,(2,j,k,mu))
           continue
         end

         mu = -Xij - Xik + Xjk
         if mu > epsi
             push!(local_corr,(3,j,k,mu))
         end
     end
     end
 end

# Given a first set of corrections that need to be made, make those corrections serially
# This may not be the most efficient way to go through the double loop, but
# either way it is so much slower than the triple loop
# it's not worth spending time optimizing it now
function FirstFix!(corrections::Vector{Vector{Tuple{Int64,Int64,Int64,Float64}}},n::Int64,W::Matrix{Float64},E::Matrix{Float64},Enew::Matrix{Float64},Elam::Matrix{Float64},lamt::Float64,Counters::Array{Int64})

    for a = 1:n
      lengtha = length(corrections[a])

      for v = 1:lengtha
          Trinum = corrections[a][Counters[a]][1]
          i = a
          j = corrections[a][v][2]
          k = corrections[a][v][3]
          mu = corrections[a][v][4]

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
          signs[Trinum] = -1;

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
end

function FirstDouble!(Fnew::Matrix{Float64},F::Matrix{Float64},Flam::Matrix{Float64},Enew::Matrix{Float64},E::Matrix{Float64},Elam::Matrix{Float64},P::Matrix{Float64},Q::Matrix{Float64},n::Int64,lamt::Float64)
    epsi = 1e-8
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

end


function MainTripleLoop!(i::Int64,D::Matrix{Int8},E::Matrix{Float64},F::Matrix{Float64},W::Matrix{Float64},Enew::Matrix{Float64},Fnew::Matrix{Float64},Elam::Matrix{Float64},Flam::Matrix{Float64},LengthI::Int64,n::Int64,L_curr::Vector{Tuple{Int64,Int64,Int64,Float64}},L_next::Vector{Tuple{Int64,Int64,Int64,Float64}})

if LengthI == 0
    return
end

println("Node $i starting the main triple loop, has list length $LengthI ")



#         if ListLength[i] == 0
#             println("Found a zero")
#             continue
#         else
#             println("Is this an issue?")
#             nextTriplet = L_curr[i][1]
#         end
#
#        for j = i+1:n-1
#        Eij = E[j,i]
#        Dij = D[j,i]
#        Wij  = W[j,i]
#         for k = j+1:n
#         Dik = D[k,i]
#         Djk = D[k,j]
#         Eik = E[k,i]
#         Ejk = E[k,j]
#         Wik = W[k,i]
#         Wjk = W[k,j]
#
#         # Consider defining a new function here that you can just pass things off
#         # to three times so that you don't have to repeat a bunch of code
#         #
#         # Need to understand how different threads will treat the fact that we're
#         # re-naming Eik, Dik, etc. every time. Will these be local to the thread
#         # that defines them?
#
#         # Check triplet i,j,k
#         if j == nextTriplet[i][2] && k == nextTriplet[i][3] && 1 == nextTriplet[i][1]
#           cor = nextTriplet[i][4]
#
#           # We need to scale the correction since we're projecting into the minimum
#           # weighted 2-norm direction
#           denom = Wij*Wjk + Wik*Wij + Wjk*Wik
#
#           Eij = Eij + cor*(Wik*Wjk/denom)
#           Eik = Eik - cor*(Wij*Wjk/denom)
#           Ejk = Ejk - cor*(Wik*Wij/denom)
#
#           # Move along in the list of triplets with corrections
#           if Counters[i] < ListLength[i]
#             Counters[i] +=1
#             nextTriplet[i] = L_curr[i][Counters[i]]
#           end
#         end
#
#         b = Dik + Djk - Dij
#         mu = (Eij - Ejk - Eik - b)
#
#         if mu > epsi
#           push!(L_next[i],(1,j,k,mu))
#           # no further adjustment will be needed for this triplet if this happened
#           continue
#         end
#         # Done checking triangle i,j,k
#
#         # Check triplet i,k,j
#         if j == nextTriplet[i][2] && k == nextTriplet[i][3] && 1 == nextTriplet[i][2]
#           cor = nextTriplet[i][4]
#
#           # We need to scale the correction since we're projecting into the minimum
#           # weighted 2-norm direction
#           denom = Wij*Wjk + Wik*Wij + Wjk*Wik
#
#           Eij = Eij - cor*(Wik*Wjk/denom)
#           Eik = Eik + cor*(Wij*Wjk/denom)
#           Ejk = Ejk - cor*(Wik*Wij/denom)
#
#           # Move along in the list of triplets with corrections
#           if Counters[i] < ListLength[i]
#             Counters[i] +=1
#             nextTriplet[i] = L_curr[i][Counters[i]]
#           end
#         end
#
#         b = -Dik + Djk + Dij
#         mu = (-Eij - Ejk + Eik - b)
#
#         if mu > epsi
#           push!(L_next[i],(2,j,k,mu))
#           # no further adjustment will be needed for this triplet if this happened
#           continue
#         end
#         # Done checking triangle i,k,j
#
#         # Check triplet j,k,i
#         if j == nextTriplet[i][2] && k == nextTriplet[i][3] && 1 == nextTriplet[i][3]
#           cor = nextTriplet[i][4]
#
#           # We need to scale the correction since we're projecting into the minimum
#           # weighted 2-norm direction
#           denom = Wij*Wjk + Wik*Wij + Wjk*Wik
#
#           Eij = Eij - cor*(Wik*Wjk/denom)
#           Eik = Eik - cor*(Wij*Wjk/denom)
#           Ejk = Ejk + cor*(Wik*Wij/denom)
#
#           # Move along in the list of triplets with corrections
#           if Counters[i] < ListLength[i]
#             Counters[i] +=1
#             nextTriplet = L_curr[i][Counters[i]]
#           end
#         end
#
#         b = Dik - Djk + Dij
#         mu = (-Eij + Ejk - Eik - b)
#
#         if mu > epsi
#           push!(L_next,(3,j,k,mu))
#           # no further adjustment will be needed for this triplet if this happened
#           continue
#         end
#         # Done checking triangle j,k,i
#
#     end
# end # end inner double for loop

end

# The main loop
function MainIteration!(D::Matrix{Int8},E::Matrix{Float64},Enew::Matrix{Float64},Elam::Matrix{Float64},F::Matrix{Float64},Fnew::Matrix{Float64},Flam::Matrix{Float64},W::Matrix{Float64},P::Matrix{Float64},Q::Matrix{Float64},L_curr::Vector{Vector{Tuple{Int64,Int64,Int64,Float64}}},L_next::Vector{Vector{Tuple{Int64,Int64,Int64,Float64}}})
println("Do the main iteration")

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

# # Initialize the nextTriplets vector, and vector of lengths of lists
# ListLength = zeros(Int64,n)
#
# # This probably won't work...need to come back later to deal with nextTriplet
# #nextTriplet = Vector{Tuple{Int64,Int64,Int64,Float64}}()
# for i = 1:n
#     #push!(nextTriplet[i],L_curr[i][1])
#     #nextTriplet[i] = L_curr[i][1]
#     ListLength[i] = length(L_curr[i])
#     ll = length(L_curr[i])
#     println("Node $i has length $ll ")
# end
#
# Counters = ones(Int64,n)


# Triple For Loop Function Here
# MainTripleLoop(params...)
# When passing through the triple loop (in parallel), we do not update the value
# of E ever, we only keep track of constraints that are violated
#
# Before checking for constraint violation, we adjust the variabled Eij, Eik,
# and Ejk if there is a nonzero correction/adjustment value in L_curr
for i = 1:n-2
    LengthI = length(L_curr[i])
    MainTripleLoop!(i,D,E,F,W,Enew,Fnew,Elam,Flam,LengthI,n,L_curr[i],L_next[i])
end


end

using MAT

function TriangleFixingStub(n)

    mat = matread("KarateA.mat")
    A = mat["A"]

    mat = matread("graphs/SmallNets.mat")
    A = mat["Afootball"]

    mat = matread("graphs/A$n\_sparse.mat")
    #A = mat["A"]

    n = size(A,1)
    Dgraph = zeros(Int8,n,n)
    for i = 1:n-1
        for j = i+1:n
            if A[i,j] < .1
                Dgraph[j,i] = 1
            end
        end
    end
    D = Dgraph

    E = zeros(n,n);
    lam = .25
    W = lam*ones(n,n)
    for i = 1:n-1
        for j = i+1:n
            if A[i,j] > .1
                W[j,i] = 1-lam
            end
        end
    end
    gam = 1;
    Etol = .05
    TriTol = .05
    F = -gam*ones(n,n)
    P = zeros(n,n)
    Q = zeros(n,n)

    @show Threads.nthreads()

    L_curr = Vector{Vector{Tuple{Int64,Int64,Int64,Float64}}}()
    for i = 1:n
        push!(L_curr,Vector{Tuple{Int64,Int64,Int64,Float64}}())
    end
    @time FirstTripleLoop!(E,D,L_curr)
    TotalViolations = 0
    for i = 1:n
        TotalViolations += length(L_curr[i])
    end
    println("Number of Initial violations: $TotalViolations")

    # STEP 2: Adjust E based on updates from triplet constraints

   # After we have collected all the things that need to be adjusted (in parallel)
   # we must go through the list and make appropriate adjustments

   # This counts through which index of corrections[i] we are at
   Counters = ones(Int64,n)

   # Initialize a zero matrix
   Enew = zeros(n,n)
   Elam = zeros(n,n)

   # each adjustment will be weighted equally among all constraints
   numcon = n*(n-1)^2/2
   lamt = 1/numcon

   # Make adjustements and add them to the new matrix G
   # This is done sequentially
   @time FirstFix!(L_curr,n,W,E,Enew,Elam,lamt,Counters)

   ## Now we handle the constraints of the form Eij - Fij <= 0,
   # -Eij - Fij <= 0.
   # We do not consider the constraints violated and projections made by
   # considering the triangle inequality consraints.
   # Here we again go back to what is currenctly stored in E and F originally.
   # We store the updates in Enew, not in E. And we scale the update by the corresponding
   # lambda weight.

   Flam = zeros(n,n)    # lambda weight for F matrix
   Fnew = zeros(n,n)    # next iterate for F matrix
   P = zeros(n,n)   # correction variables for constraints -Eij-Fij <=0
   Q = zeros(n,n)  # correction variables for constraints Eij-Fij <=0

   @time FirstDouble!(Fnew,F,Flam,Enew,E,Elam,P,Q,n,lamt)

   # In future rounds, L_current stores corrections that need to be made before
   # making projections, while L_next keeps track of new constraint violations
   # we need to keep track of
   L_next = Vector{Vector{Tuple{Int64,Int64,Int64,Float64}}}()
   for i = 1:n
       push!(L_next,Vector{Tuple{Int64,Int64,Int64,Float64}}())
   end

   # Now we enter the main loop

   maxits = 1

   for iter = 1:maxits
       tempE = E[:,:]

       # We are ready to reset E and Enew, F and Fnew
       F = Fnew
       Fnew = zeros(n,n)
       Flam = zeros(n,n)

       E = Enew
       Enew = zeros(n,n)
       Elam = zeros(n,n)
       MainIteration!(D,E,Enew,Elam,F,Fnew,Flam,W,P,Q,L_curr,L_next)

       # Get rid of old corrections, prepare space for new
       L_curr = L_next
       L_next = Vector{Vector{Tuple{Int64,Int64,Int64,Float64}}}()
       for i = 1:n
           push!(L_next,Vector{Tuple{Int64,Int64,Int64,Float64}}())
       end

   end  # end main loo

end # end TriangleFixingStub

TriangleFixingStub(parse(Int64,ARGS[1]))
