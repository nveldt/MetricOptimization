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

# Go through the triple loop the first time, and identify
# triangle inequality violations that must be fixed
function FirstTripleLoop!(E::Matrix{Float64},D::Matrix{Int8},corrections::Vector{Vector{Tuple{Int64,Int64,Int64,Float64,Float64}}})
    n = size(E,1)
     Threads.@threads for j = 2:n-1
        JobJ!(E,D,j,corrections[j])
    end
end

function JobJ!(E::Matrix{Float64},D::Matrix{Int8},j::Int64,local_corr::Vector{Tuple{Int64,Int64,Int64,Float64,Float64}})
    # This is a task performed by a single thread. It is the inner two loops
    # of the triple for loop
    n = size(E,1); epsi = 1e-8
    for i = 1:j-1
         Xij = E[j,i] + D[j,i]
         for k = j+1:n
         Xik = E[k,i] + D[k,i]; Xjk = E[k,j] + D[k,j]

         mu = Xij - Xik - Xjk
         if mu > epsi
           push!(local_corr,(1,i,k,0,mu))
           continue
         end

         mu = -Xij + Xik - Xjk
         if mu > epsi
           push!(local_corr,(2,i,k,0,mu))
           continue
         end

         mu = -Xij - Xik + Xjk
         if mu > epsi
             push!(local_corr,(3,i,k,0,mu))
         end
     end
     end
 end

# Given a first set of corrections that need to be made, make those corrections serially
# This may not be the most efficient way to go through the double loop, but
# either way it is so much slower than the triple loop
# it's not worth spending time optimizing it now
function FirstFix!(corrections::Vector{Vector{Tuple{Int64,Int64,Int64,Float64,Float64}}},n::Int64,W::Matrix{Float64},E::Matrix{Float64},Enew::Matrix{Float64},Elam::Matrix{Float64},lamt::Float64)

    # corrections is now indexed by j, the second largest index in the triplet i, j, k
    # I.e. corrections[j] stores all triplets of the form (i,j,k) where i<j<k
    # such that one of the three triangle inequality constraints at this triplet
    # was not satisfied in a recent pass through the constraints.
    for a = 2:n-1
      lengtha = length(corrections[a])

      for v = 1:lengtha
          Trinum = corrections[a][v][1]
          j = a
          i = corrections[a][v][2]
          k = corrections[a][v][3]
          mu = corrections[a][v][5]

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

          # Now since there was an adjustment, the "adjustment" gets moved
          # to the "correction" part of the tuple
          corrections[a][v] = (Trinum,i,k,mu,0.0)
      end
    end
end

function MainTripleLoop!(j::Int64,D::Matrix{Int8},E::Matrix{Float64},W::Matrix{Float64},LengthJ::Int64,n::Int64,l_curr::Vector{Tuple{Int64,Int64,Int64,Float64,Float64}},l_next::Vector{Tuple{Int64,Int64,Int64,Float64,Float64}},updates::Vector{Vector{Tuple{Int64,Int64,Int64,Float64}}})
    epsi = 1e-8

    # Go through loop,
    #   - make corrections based on nonzero tuples in l_curr (= L_curr[j])
    #   - store new violation information in l_next (= L_next[j])
    #   - if a contraint has a nonzero correction, but then no violation,
    #       this still would affect the update of E eventually, so we store it
    #       in the 'updates' vector

    # Keep track of where we are in l_curr
    TripInd = 1
    nextTriplet = l_curr[TripInd]

    # nextTriplet has the form (TriNum, i, k, cor, 0)
    # indicating that triangle (i,j,k) needs a correction 'cor' at constraint TriNum

    for i = 1:j-1
        Eij = E[j,i]
        Dij = D[j,i]
        Wij = W[j,i]

        for k = j+1:n
            Eik = E[k,i]
            Dik = D[k,i]
            Ejk = E[k,j]
            Djk = D[k,j]
            Wik = W[k,i]
            Wjk = W[k,j]

      if i == nextTriplet[2] && k == nextTriplet[3] # This means one of three constraints corresponds to a correction

            # Check constraint (1): Xij - Xik - Xjk <= 0, (recall Eab = Xab - Dab)
            denom = Wij*Wjk + Wik*Wij + Wjk*Wik

            if i == nextTriplet[2] && k == nextTriplet[3] && 1 == nextTriplet[1]
              # Scale to project with minimum W-norm distance
              cor = nextTriplet[4]
              Eij = Eij + cor*(Wik*Wjk/denom)
              Eik = Eik - cor*(Wij*Wjk/denom)
              Ejk = Ejk - cor*(Wik*Wij/denom)

              # Move along in the list of triplets with corrections
              if TripInd < LengthJ
                TripInd +=1
                nextTriplet = l_curr[TripInd]
              end
            else
              cor = 0
            end

            b = Dik + Djk - Dij
            mu = (Eij - Ejk - Eik - b)
            if mu > epsi
              push!(l_next,(1,i,k,cor,mu))
              continue
            end
            if mu <= epsi && cor > epsi
                # Then we have an update at this constraint but not a future correction
                push!(updates,(1,i,k,cor))
            end
            # Done with Xij - Xik - Xjk <= 0

            #Check constraint (2): -Xij + Xik - Xjk <= 0
            if i == nextTriplet[2] && k == nextTriplet[3] && 2 == nextTriplet[1]
              cor = nextTriplet[4]
              # Scale to project with minimum W-norm distance
              Eij = E[j,i] - cor*(Wik*Wjk/denom)
              Eik = E[k,i] + cor*(Wij*Wjk/denom)
              Ejk = E[k,j] - cor*(Wik*Wij/denom)

              # Move along in the list of triplets with corrections
              if TripInd < LengthJ
                TripInd +=1
                nextTriplet = l_curr[TripInd]
              end
            else
                Eij = E[j,i]
                Eik = E[k,i]
                Ejk = E[k,j]
                cor = 0
            end

            b = -Dik + Djk + Dij
            mu = (-Eij - Ejk + Eik - b)
            if mu > epsi
              push!(l_next,(2,i,k,cor,mu))
              continue
            end
            if mu <= epsi && cor > epsi
                # Then we have an update at this constraint but not a future correction
                push!(updates,(2,i,k,cor))
            end
            # Done checking (2): -Xij + Xik - Xjk <= 0

            #Check constraint (3): -Xij - Xik + Xjk <= 0
            if i == nextTriplet[2] && k == nextTriplet[3] && 3 == nextTriplet[1]  #nonzero correction
              cor = nextTriplet[4]
              # Scale to project with minimum W-norm distance
              Eij = E[j,i] - cor*(Wik*Wjk/denom)
              Eik = E[k,i] - cor*(Wij*Wjk/denom)
              Ejk = E[k,j] + cor*(Wik*Wij/denom)

              # Move along in the list of triplets with corrections
              if TripInd < LengthJ
                TripInd +=1
                nextTriplet = l_curr[TripInd]
              end
            else                                                                  # zero correction
                Eij = E[j,i]
                Eik = E[k,i]
                Ejk = E[k,j]
                cor = 0
            end

            b = Dik - Djk + Dij
            mu = (-Eij + Ejk - Eik - b)
            if mu > epsi                    # non-trivial projection
              push!(l_next,(3,i,k,cor,mu))
            end
            if mu <= epsi && cor > epsi     # trivial projection, nontrivial correction.
              push!(updates,(3,i,k,cor))
            end


      else        # This is if we don't have any corrections for this triplet at all
          Xij = E[j,i] + D[j,i]
          Xik = E[k,i] + D[k,i]; Xjk = E[k,j] + D[k,j]

         mu = Xij - Xik - Xjk
         if mu > epsi
           push!(l_next,(1,i,k,0.0,mu))
           continue
         end

         mu = -Xij + Xik - Xjk
         if mu > epsi
           push!(l_next,(2,i,k,0.0,mu))
           continue
         end

         mu = -Xij - Xik + Xjk
         if mu > epsi
             push!(l_next,(3,i,k,0.0,mu))
         end
     end

    end # end inner double loop 1
    end # end inner double loop 2
end # End MainTripleLoop

function Double_Loop!(Fnew::Matrix{Float64},F::Matrix{Float64},Flam::Matrix{Float64},Enew::Matrix{Float64},E::Matrix{Float64},Elam::Matrix{Float64},P::Matrix{Float64},Q::Matrix{Float64},n::Int64,lamt::Float64)
    epsi = 1e-8
    for i = 1:n-1
      for j = i+1:n

        cor = P[j,i]
        Eij = E[j,i] - cor
        Fij = F[j,i] - cor
        delta = -Eij - Fij
        if delta > epsi
          Enew[j,i] += lamt*(Eij + delta/2)
          Fnew[j,i] += lamt*(Fij + delta/2)
          P[j,i] = delta/2
        else
          P[j,i] = 0.0
        end
    end
    end

    for i = 1:n-1
      for j = i+1:n

        cor = Q[j,i]
        Eij = E[j,i] - cor
        Fij = F[j,i] + cor
        delta = Eij - Fij
        if delta > epsi
          Enew[j,i] += lamt*(Eij - delta/2)
          Fnew[j,i] += lamt*(Fij + delta/2)
          Q[j,i] = -delta/2
        else
            Q[j,i] = 0.0
        end
      end
    end

    # Up until now we have only incremented entry (i,j) of Enew and Fnew
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

function UpdateE_TriConstraints!(updates::Vector{Vector{Tuple{Int64,Int64,Int64,Float64}}},L_next::Vector{Vector{Tuple{Int64,Int64,Int64,Float64,Float64}}},n::Int64,W::Matrix{Float64},E::Matrix{Float64},Enew::Matrix{Float64},Elam::Matrix{Float64},lamt::Float64)
    # L_next[j] stores tuples of the form (Trinum,i,k,cor,mu)
    #   where Trinum, i,k identifies a constraint with a nonzero mu (and possibly
    # nonzero cor), and the idea is we need to update things
    for a = 2:n-1
      lengtha = length(L_next[a])

      for v = 1:lengtha
          Trinum = L_next[a][v][1]
          j = a
          i = L_next[a][v][2]
          k = L_next[a][v][3]
          cor = L_next[a][v][4]
          mu = L_next[a][v][5]

          # This means there's a triplet (i,j,k), and Trinum tells you which variable
          # Eij, Eik, Ejk is the "odd one out" in the constraint.
          Wik = W[k,i]
          Wjk = W[k,j]
          Wij = W[j,i]
          Eik = E[k,i]
          Ejk = E[k,j]
          Eij = E[j,i]

          denom = Wij*Wjk + Wik*Wij + Wjk*Wik

          signs = -ones(3,1)
          signs[Trinum] = 1;
          # Two of them have an added term, and the other has a term subtracted
          Enew[j,i] += lamt*(Eij + signs[1]*(cor-mu)*(Wik*Wjk/denom))
          Enew[k,i] += lamt*(Eik + signs[2]*(cor-mu)*(Wij*Wjk/denom))
          Enew[k,j] += lamt*(Ejk + signs[3]*(cor-mu)*(Wik*Wij/denom))

          # Update how much "lambda-weight" has been added to each entry
          Elam[j,i] += lamt
          Elam[k,i] += lamt
          Elam[k,j] += lamt

          # Now since there was an adjustment, the "adjustment" gets moved
          # to the "correction" part of the tuple
          L_next[a][v] = (Trinum,i,k,mu,0)
      end
    end


end # End UpdateE_TriConstraints

# Skeleton for the main loop
function MainIteration!(D::Matrix{Int8},E::Matrix{Float64},Enew::Matrix{Float64},Elam::Matrix{Float64},F::Matrix{Float64},Fnew::Matrix{Float64},Flam::Matrix{Float64},W::Matrix{Float64},P::Matrix{Float64},Q::Matrix{Float64},L_curr::Vector{Vector{Tuple{Int64,Int64,Int64,Float64,Float64}}},L_next::Vector{Vector{Tuple{Int64,Int64,Int64,Float64,Float64}}},updates::Vector{Vector{Tuple{Int64,Int64,Int64,Float64}}})
    #println("Do the main iteration")
    n = size(D,1)
    epsi = 1e-8                 # constraint tolerance
    numcon = n*(n-1)^2/2        # number of constraints
    lamt = 1/numcon

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

    # STEP 1: iterate through all triplets and keep track of
    #         nontrivial corrections and nontrivial projections
  for j = 2:n-1
        LengthJ = length(L_curr[j]) # total number of corrections thread j must perform
        if LengthJ == 0
            # Have it taken care of by a loop that doesn't require adjustments
            JobJ!(E,D,j,L_next[j])
        else
            MainTripleLoop!(j,D,E,W,LengthJ,n,L_curr[j],L_next[j],updates)
        end
    end

    # STEP 2: Go through L_next, make appropriate adjustments to E, and
    #         prepare L_next to be the new L_curr
    #         Then visit each thing in 'updates' and perform updates based on that too
    UpdateE_TriConstraints!(updates,L_next,n,W,E,Enew,Elam,lamt)

    # STEP 3: Visit the double loop with constraints of the form
    #           eij - fij <= 0 and -eij - fij <= 0
    #           and perform correction, projections, and updates all at once
    Double_Loop!(Fnew,F,Flam,Enew,E,Elam,P,Q,n,lamt)

end

using MAT

function TriangleFixingStub(n)

    mat = matread("KarateA.mat")
    A = mat["A"]

    # mat = matread("graphs/SmallNets.mat")
    # A = mat["Afootball"]
    #
    # mat = matread("graphs/A$n\_sparse.mat")
    # A = mat["A"]

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
    lam = .9
    W = lam*ones(n,n)
    for i = 1:n-1
        for j = i+1:n
            if A[i,j] > .1
                W[j,i] = 1-lam
            end
        end
    end
    gam = 5
    Etol = .05
    TriTol = .01
    F = -gam*ones(n,n)
    P = zeros(n,n)
    Q = zeros(n,n)

    @show Threads.nthreads()

    L_curr = Vector{Vector{Tuple{Int64,Int64,Int64,Float64,Float64}}}()
    for i = 1:n
        push!(L_curr,Vector{Tuple{Int64,Int64,Int64,Float64,Float64}}())
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
   @time FirstFix!(L_curr,n,W,E,Enew,Elam,lamt)

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

   @time Double_Loop!(Fnew,F,Flam,Enew,E,Elam,P,Q,n,lamt)

   # In future rounds, L_current stores corrections that need to be made before
   # making projections, while L_next keeps track of new constraint violations
   # we need to keep track of
   L_next = Vector{Vector{Tuple{Int64,Int64,Int64,Float64,Float64}}}()
   for j = 1:n
       push!(L_next,Vector{Tuple{Int64,Int64,Int64,Float64,Float64}}())
   end

   # updates keeps track of constraints where a nontrivial correction terms exists
   # but then there is no nontrivial projection
   updates = Vector{Vector{Tuple{Int64,Int64,Int64,Float64}}}()
   for j = 1:n
       push!(updates,Vector{Tuple{Int64,Int64,Int64,Float64}}())
   end

   # Now we enter the main loop

   # We are ready to reset E and Enew, F and Fnew
   F = Fnew
   Fnew = zeros(n,n)
   Flam = zeros(n,n)

   E = Enew
   Enew = zeros(n,n)
   Elam = zeros(n,n)

iter = 0
   while true
     iter +=1
       tempE = E[:,:]
       tempF = F[:,:]

      tic()
      MainIteration!(D,E,Enew,Elam,F,Fnew,Flam,W,P,Q,L_curr,L_next,updates)
      lasttime = toq()
       TotalViolations = 0
       for j = 1:n
           TotalViolations += length(L_next[j])
       end

       # Get rid of old corrections, prepare space for new
       L_curr = L_next

       L_next = Vector{Vector{Tuple{Int64,Int64,Int64,Float64,Float64}}}()
       for j = 1:n
           push!(L_next,Vector{Tuple{Int64,Int64,Int64,Float64,Float64}}())
       end

      updates = Vector{Vector{Tuple{Int64,Int64,Int64,Float64}}}()
      for j = 1:n
          push!(updates,Vector{Tuple{Int64,Int64,Int64,Float64}}())
      end
       # We are ready to reset E and Enew, F and Fnew
       F = Fnew
       Fnew = zeros(n,n)
       Flam = zeros(n,n)

       E = Enew
       Enew = zeros(n,n)
       Elam = zeros(n,n)

       change = norm(vec(E-tempE)) + norm(vec(F-tempF))
       tricheck = TriangleCheck(E+D,TriTol)
       objective = LPcc_obj(A,E+D,lam)
       ch = round(change,3)
       lt = round(lasttime,2)
       ob = round(objective,3)
       if iter%10 == 1
         G = E+D
         tr = FullTriangleCheck(G)
         println("Triangle inequality check: $tr")
       end
       println("Iteration $iter, E changed by $ch, Numvio = $TotalViolations 3Loop: $lt, Obj: $ob")
       if change < Etol && tricheck
         break
       end
   end  # end main loop
   #@show E+D
end # end TriangleFixingStub

TriangleFixingStub(parse(Int64,ARGS[1]))
