# Go through the triple loop the first time, and identify
# triangle inequality violations that must be fixed
function FirstTripleLoop!(E::Matrix{Float64},D::Matrix{Int8},corrections::Vector{Vector{Tuple{Int64,Int64,Int64,Float64}}})
    n = size(E,1)
     Threads.@threads for j = 2:n-1
        JobJ!(E,D,j,corrections[j])
    end
end

function JobJ!(E::Matrix{Float64},D::Matrix{Int8},j::Int64,local_corr::Vector{Tuple{Int64,Int64,Int64,Float64}})
    # This is a task performed by a single thread. It is the inner two loops
    # of the triple for loop
    n = size(E,1); epsi = 1e-8
    for i = 1:j-1
        Xij = E[j,i] + D[j,i]
         for k = j+1:n
         Xik = E[k,i] + D[k,i]; Xjk = E[k,j] + D[k,j]

         mu = Xij - Xik - Xjk
         if mu > epsi
           push!(local_corr,(1,i,k,mu))
           continue
         end

         mu = -Xij + Xik - Xjk
         if mu > epsi
           push!(local_corr,(2,i,k,mu))
           continue
         end

         mu = -Xij - Xik + Xjk
         if mu > epsi
             push!(local_corr,(3,i,k,mu))
         end
     end
     end
 end

# Given a first set of corrections that need to be made, make those corrections serially
# This may not be the most efficient way to go through the double loop, but
# either way it is so much slower than the triple loop
# it's not worth spending time optimizing it now
function FirstFix!(corrections::Vector{Vector{Tuple{Int64,Int64,Int64,Float64}}},n::Int64,W::Matrix{Float64},E::Matrix{Float64},Enew::Matrix{Float64},Elam::Matrix{Float64},lamt::Float64,Counters::Array{Int64})

    # corrections is now indexed by j, the second largest index in the triplet i, j, k
    # I.e. corrections[j] stores all triplets of the form (i,j,k) where i<j<k
    # such that one of the three triangle inequality constraints at this triplet
    # was not satisfied in a recent pass through the constraints.
    for a = 2:n-1
      lengtha = length(corrections[a])

      for v = 1:lengtha
          Trinum = corrections[a][Counters[a]][1]
          j = a
          i = corrections[a][v][2]
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

function MainTripleLoop2!(j::Int64,D::Matrix{Int8},E::Matrix{Float64},F::Matrix{Float64},W::Matrix{Float64},Enew::Matrix{Float64},Fnew::Matrix{Float64},Elam::Matrix{Float64},Flam::Matrix{Float64},LengthJ::Int64,n::Int64,l_curr::Vector{Tuple{Int64,Int64,Int64,Float64}},l_next::Vector{Tuple{Int64,Int64,Int64,Float64}})
    epsi = 1e-8

    # Go through loop,
    #   - make corrections based on nonzero tuples in l_curr (= L_curr[j])
    #   - store new violation information in l_next (= L_next[j])

    # Keep track of where we are in l_curr
    # nextTriplet = l_curr[1]

    nT1 = l_curr[1][1]
    nT2 = l_curr[1][2]
    nT3 = l_curr[1][3]
    cor = l_curr[1][4]

    TripInd = 1

    # nextTriplet has the form (TriNum, i, k, cor)
    # indicating that triangle (i,j,k) needs a correction 'cor' at constraint TriNum

    for i = 1:j-1
        Eij = E[j,i]
        Dij = D[j,i]
        Wij  = W[j,i]
        Oij = Eij

        for k = j+1:n
            Eik = E[k,i]
            Dik = D[k,i]
            Ejk = E[k,j]
            Djk = D[k,j]

            Oik = Eik
            Ojk = Ejk

        if i == nT2 && k == nT3
            # Check constraint (1): Xij - Xik - Xjk <= 0, (recall Eab = Xab - Dab)
            #if 1 == nextTriplet[1] && i == nextTriplet[2] && k == nextTriplet[3]
            Wik = W[k,i]; Wjk = W[k,j]
            denom = Wij*Wjk + Wik*Wij + Wjk*Wik

            if nT1 == 1
              # Scale to project with minimum W-norm distance
              Eij = Oij + cor*(Wik*Wjk/denom)
              Eik = Oik - cor*(Wij*Wjk/denom)
              Ejk = Ojk - cor*(Wik*Wij/denom)

              # Move along in the list of triplets with corrections
              if TripInd < LengthJ
                #nextTriplet = l_curr[TripInd]
                TripInd +=1
                nT1 = l_curr[TripInd][1]
                nT2 = l_curr[TripInd][2]
                nT3 = l_curr[TripInd][3]
                cor = l_curr[TripInd][4]
              end
          end

            b = Dik + Djk - Dij
            mu = (Eij - Ejk - Eik - b)
            if mu > epsi
              push!(l_next,(1,i,k,mu))
              continue
            end
            # Done with Xij - Xik - Xjk <= 0

            #Check constraint (2): -Xij + Xik - Xjk <= 0
            if nT1 == 2
              Eij = Oij - cor*(Wik*Wjk/denom)
              Eik = Oik + cor*(Wij*Wjk/denom)
              Ejk = Ojk - cor*(Wik*Wij/denom)

              # Move along in the list of triplets with corrections
              if TripInd < LengthJ
                TripInd +=1
                nT1 = l_curr[TripInd][1]
                nT2 = l_curr[TripInd][2]
                nT3 = l_curr[TripInd][3]
                cor = l_curr[TripInd][4]
              end
            else
                Eij = Oij
                Eik = Oik
                Ejk = Ojk
            end

            b = -Dik + Djk + Dij
            mu = (-Eij - Ejk + Eik - b)
            if mu > epsi
              push!(l_next,(2,i,k,mu))
              continue
            end
            # Done checking (2): -Xij + Xik - Xjk <= 0

            #Check constraint (3): -Xij - Xik + Xjk <= 0
            if nT1 == 3
              Eij = Oij - cor*(Wik*Wjk/denom)
              Eik = Oik - cor*(Wij*Wjk/denom)
              Ejk = Ojk + cor*(Wik*Wij/denom)

              # Move along in the list of triplets with corrections
              if TripInd < LengthJ
                TripInd +=1
                nT1 = l_curr[TripInd][1]
                nT2 = l_curr[TripInd][2]
                nT3 = l_curr[TripInd][3]
                cor = l_curr[TripInd][4]
              end
            else
              Eij = Oij
              Eik = Oik
              Ejk = Ojk
            end

            b = Dik - Djk + Dij
            mu = (-Eij + Ejk - Eik - b)
            if mu > epsi
              push!(l_next,(3,i,k,mu))
            end

    else
        Xij = E[j,i] + D[j,i]
        Xik = E[k,i] + D[k,i]; Xjk = E[k,j] + D[k,j]

         mu = Xij - Xik - Xjk
         if mu > epsi
           push!(l_next,(1,i,k,mu))
           continue
         end

         mu = -Xij + Xik - Xjk
         if mu > epsi
           push!(l_next,(2,i,k,mu))
           continue
         end

         mu = -Xij - Xik + Xjk
         if mu > epsi
             push!(l_next,(3,i,k,mu))
         end
     end

    end # end inner double loop 1
    end # end inner double loop 2
end # End MainTripleLoop2

function MainTripleLoop!(j::Int64,D::Matrix{Int8},E::Matrix{Float64},F::Matrix{Float64},W::Matrix{Float64},Enew::Matrix{Float64},Fnew::Matrix{Float64},Elam::Matrix{Float64},Flam::Matrix{Float64},LengthJ::Int64,n::Int64,l_curr::Vector{Tuple{Int64,Int64,Int64,Float64}},l_next::Vector{Tuple{Int64,Int64,Int64,Float64}})
    epsi = 1e-8
    #println("Node $j starting the main triple loop, has list length $LengthJ ")

    # Go through loop,
    #   - make corrections based on nonzero tuples in l_curr (= L_curr[j])
    #   - store new violation information in l_next (= L_next[j])

    # Keep track of where we are in l_curr
    nextTriplet = l_curr[1]
    TripInd = 1

    # nextTriplet has the form (TriNum, i, k, cor)
    # indicating that triangle (i,j,k) needs a correction 'cor' at constraint TriNum

    for i = 1:j-1
        Eij = E[j,i]
        Dij = D[j,i]
        Wij  = W[j,i]
        Oij = Eij

        for k = j+1:n
            Eik = E[k,i]
            Dik = D[k,i]
            Ejk = E[k,j]
            Djk = D[k,j]
            Wik = W[k,i]
            Wjk = W[k,j]

            Oik = Eik
            Ojk = Ejk

            # Check constraint (1): Xij - Xik - Xjk <= 0, (recall Eab = Xab - Dab)
            if 1 == nextTriplet[1] && i == nextTriplet[2] && k == nextTriplet[3]
              cor = nextTriplet[4]

              # Scale to project with minimum W-norm distance
              denom = Wij*Wjk + Wik*Wij + Wjk*Wik

              Eij = Oij + cor*(Wik*Wjk/denom)
              Eik = Oik - cor*(Wij*Wjk/denom)
              Ejk = Ojk - cor*(Wik*Wij/denom)

              # Move along in the list of triplets with corrections
              if TripInd < LengthJ
                TripInd +=1
                nextTriplet = l_curr[TripInd]
              end
          end

            b = Dik + Djk - Dij
            mu = (Eij - Ejk - Eik - b)
            if mu > epsi
              push!(l_next,(1,i,k,mu))
              continue
            end
            # Done with Xij - Xik - Xjk <= 0

            #Check constraint (2): -Xij + Xik - Xjk <= 0
            if i == nextTriplet[2] && k == nextTriplet[3] && 2 == nextTriplet[1]
              cor = nextTriplet[4]

              # Scale to project with minimum W-norm distance
              denom = Wij*Wjk + Wik*Wij + Wjk*Wik

              Eij = Oij - cor*(Wik*Wjk/denom)
              Eik = Oik + cor*(Wij*Wjk/denom)
              Ejk = Ojk - cor*(Wik*Wij/denom)

              # Move along in the list of triplets with corrections
              if TripInd < LengthJ
                TripInd +=1
                nextTriplet = l_curr[TripInd]
              end
            else
                Eij = Oij
                Eik = Oik
                Ejk = Ojk
            end

            b = -Dik + Djk + Dij
            mu = (-Eij - Ejk + Eik - b)
            if mu > epsi
              push!(l_next,(2,i,k,mu))
              continue
            end
            # Done checking (2): -Xij + Xik - Xjk <= 0

            #Check constraint (3): -Xij - Xik + Xjk <= 0
            if i == nextTriplet[2] && k == nextTriplet[3] && 3 == nextTriplet[1]
              cor = nextTriplet[4]

              # Scale to project with minimum W-norm distance
              denom = Wij*Wjk + Wik*Wij + Wjk*Wik

              Eij = Oij - cor*(Wik*Wjk/denom)
              Eik = Oik - cor*(Wij*Wjk/denom)
              Ejk = Ojk + cor*(Wik*Wij/denom)

              # Move along in the list of triplets with corrections
              if TripInd < LengthJ
                TripInd +=1
                nextTriplet = l_curr[TripInd]
              end
            else
              Eij = Oij
              Eik = Oik
              Ejk = Ojk
            end

            b = Dik - Djk + Dij
            mu = (-Eij + Ejk - Eik - b)
            if mu > epsi
              push!(l_next,(3,i,k,mu))
            end

            end
        end # end inner double loop

end # End MainTripleLoop

# Skeleton for the main loop
function MainIteration!(D::Matrix{Int8},E::Matrix{Float64},Enew::Matrix{Float64},Elam::Matrix{Float64},F::Matrix{Float64},Fnew::Matrix{Float64},Flam::Matrix{Float64},W::Matrix{Float64},P::Matrix{Float64},Q::Matrix{Float64},L_curr::Vector{Vector{Tuple{Int64,Int64,Int64,Float64}}},L_next::Vector{Vector{Tuple{Int64,Int64,Int64,Float64}}})
    println("Do the main iteration")
    n = size(D,1)
    epsi = 1e-8                 # constraint tolerance
    numcon = n*(n-1)^2/2        # number of constraints

    # L_curr[j] stores nonzero adjustments needed for triplets (i,j,k) with i < j < k
    # L_next[j] is an empty vector of tuples that will store new adjustments needed to violated constraints

    # L_next[j] is of tuples of the form (Int, Int, Int, Float)
    # e.g. (1,i,k, mu)
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
    Threads.@threads for j = 2:n-1
        LengthJ = length(L_curr[j])
        if LengthJ == 0
            # Have it taken care of by a loop that doesn't require adjustments
            JobJ!(E,D,j,L_next[j])
        else
            MainTripleLoop2!(j,D,E,F,W,Enew,Fnew,Elam,Flam,LengthJ,n,L_curr[j],L_next[j])
        end

    end
end



using MAT

function TriangleFixingStub(n)

    mat = matread("KarateA.mat")
    A = mat["A"]

    mat = matread("graphs/SmallNets.mat")
    A = mat["Afootball"]

    mat = matread("graphs/A$n\_sparse.mat")
    A = mat["A"]

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
    gam = 1
    Etol = .5
    TriTol = .5
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
   for j = 1:n
       push!(L_next,Vector{Tuple{Int64,Int64,Int64,Float64}}())
   end

   # Now we enter the main loop

   # We are ready to reset E and Enew, F and Fnew
   F = Fnew
   Fnew = zeros(n,n)
   Flam = zeros(n,n)

   E = Enew
   Enew = zeros(n,n)
   Elam = zeros(n,n)

   maxits = 1

   for iter = 1:maxits
       tempE = E[:,:]

       TotalViolations = 0
       for j = 1:n
           TotalViolations += length(L_next[j])
       end
       println("Before main iteration: $TotalViolations")

       @time MainIteration!(D,E,Enew,Elam,F,Fnew,Flam,W,P,Q,L_curr,L_next)

       TotalViolations = 0
       for j = 1:n
           TotalViolations += length(L_next[j])
       end
       println("Number of new violations: $TotalViolations")

       # Get rid of old corrections, prepare space for new
       L_curr = L_next
       L_next = Vector{Vector{Tuple{Int64,Int64,Int64,Float64}}}()
       for j = 1:n
           push!(L_next,Vector{Tuple{Int64,Int64,Int64,Float64}}())
       end
       # We are ready to reset E and Enew, F and Fnew
       F = Fnew
       Fnew = zeros(n,n)
       Flam = zeros(n,n)

       E = Enew
       Enew = zeros(n,n)
       Elam = zeros(n,n)



   end  # end main loo

end # end TriangleFixingStub

TriangleFixingStub(parse(Int64,ARGS[1]))
