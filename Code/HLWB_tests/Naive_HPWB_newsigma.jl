# This is a naive implementation of the Halpern-Lions_Wittmann_Buschke algorithm
# just to see if it might converge

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
          write(f, "Writing output from calling the Triangle Fixing Algorithm, HLWB naive style\n")
          write(f, "Lambda = $lam, gamma = $gam, Etol = $Etol, TriTol = $TriTol \n")
  end

# Tolerance for inequality constraints
epsi = 0
tr = 0

Kcounter = 0;

# Part 1: Triangle inequality constraint checking
@inbounds for i = 1:n-2
   for j = i+1:n-1
   Eij = E[j,i]
   Dij = D[j,i]
   Wij = W[j,i]
    for k = j+1:n

    Kcounter += 1;
    sigmaK = 1/((Kcounter + 1))
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

      continue
    end
    # Done checking triangle j,k,i


    # Now do the HLWB step with the steering sequence
    E = (1-sigmaK)*E
    F = -gam*ones(n,n)*sigmaK + (1-sigmaK)*F

end
end
end # end triple for loop


# Part 2: -E <= F and E <= F checking, the second set of constraints
for i = 1:n-1
  for j = i+1:n
    Kcounter += 1;
    Eij = E[j,i]
    Fij = F[j,i]
    delta = -Eij - Fij
    if delta > epsi
      E[j,i] = Eij + delta/2
      F[j,i] = Fij + delta/2
      P[j,i] = delta/2
      # else the correction next time will just be P[i,j] == 0
    end

    # Now do the HLWB step with the steering sequence
    sigmaK = 1/((Kcounter + 1))
    E = (1-sigmaK)*E
    F = -gam*ones(n,n)*sigmaK + (1-sigmaK)*F

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

    # Now do the HLWB step with the steering sequence
    Kcounter += 1;
    sigmaK = 1/((Kcounter + 1))
    E = (1-sigmaK)*E
    F = -gam*ones(n,n)*sigmaK + (1-sigmaK)*F


  end
end
# Done with the first time through the constraints

iter = 0
maxits = 3000
while true && iter < maxits
    tempE = E[:,:]

    # Start again with the triangle inequality constraints
    tic()
    @inbounds for i = 1:n-2
       for j = i+1:n-1
       Eij = E[j,i]
       Dij = D[j,i]
       Wij = W[j,i]
        for k = j+1:n

        Kcounter += 1;
        sigmaK = 1/((Kcounter + 1))
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

          continue
        end
        # Done checking triangle j,k,i


        # Now do the HLWB step with the steering sequence
        E = (1-sigmaK)*E
        F = -gam*ones(n,n)*sigmaK + (1-sigmaK)*F

    end
    end
    end # end triple for loop
    lasttime = toq()


    ### Now we check the non-triangle constraints
    # Part 2: -E <= F and E <= F checking

        # Part 2: -E <= F and E <= F checking, the second set of constraints
        for i = 1:n-1
          for j = i+1:n
            Kcounter += 1;
            Eij = E[j,i]
            Fij = F[j,i]
            delta = -Eij - Fij
            if delta > epsi
              E[j,i] = Eij + delta/2
              F[j,i] = Fij + delta/2
              P[j,i] = delta/2
              # else the correction next time will just be P[i,j] == 0
            end

            # Now do the HLWB step with the steering sequence
            sigmaK = 1/((Kcounter + 1))
            E = (1-sigmaK)*E
            F = -gam*ones(n,n)*sigmaK + (1-sigmaK)*F

        end
        end

        for i = 1:n-1
          for j = i+1:n
            Kcounter += 1;
            Eij = E[j,i]
            Fij = F[j,i]
            delta = Eij - Fij
            if delta > epsi
              E[j,i] = Eij - delta/2
              F[j,i] = Fij + delta/2
              Q[j,i] = -delta/2
              # else the correction next time will just be P[i,j] == 0
            end

            # Now do the HLWB step with the steering sequence
            sigmaK = 1/((Kcounter + 1))
            E = (1-sigmaK)*E
            F = -gam*ones(n,n)*sigmaK + (1-sigmaK)*F

          end
        end


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
    println("Iteration $iter, E changed by $ch, 3Loop: $lt, Obj: $ob, Last TriCheck = $tr")

    open(filename, "a") do f
      write(f,"Iteration $iter, E changed by $ch, 3Loop: $lt, Obj: $ob, Last TriCheck = $tr\n")
    end

    if change < Etol && tricheck
      break
    end

end #end while loop


return E+D


end
