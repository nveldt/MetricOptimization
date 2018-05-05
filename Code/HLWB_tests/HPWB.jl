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
# This solves the LambdaCC LP relaxation using a Triangle Fixing Algorithm
# that operates with the HLWB algorithm underneath instead of Dykstra's projection algorithm
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
epsi = 1e-12
tr = 0

Kcounter = 0;

# These keep track of the last time we updated Eij and Fij
CE = zeros(n,n);
CF = zeros(n,n);

# Part 1: Triangle inequality constraint checking
@inbounds for i = 1:n-2
   for j = i+1:n-1
    for k = j+1:n

    Kcounter += 1
    sigmaK = 1/(Kcounter + 1)

    # Update Eij, Eik, Ejk for all the updates they missed of the form Eij = (1-sigmaK)*Eij
    E[j,i] = (CE[j,i]+1)/Kcounter*E[j,i]
    E[k,i] = (CE[k,i]+1)/Kcounter*E[k,i]
    E[k,j] = (CE[k,j]+1)/Kcounter*E[k,j]

    CE[j,i] = Kcounter
    CE[k,i] = Kcounter
    CE[k,j] = Kcounter

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


    # Now do the HLWB step with the steering sequence
    E[j,i] = (Kcounter)/(Kcounter+1)*E[j,i]
    E[k,i] = (Kcounter)/(Kcounter+1)*E[k,i]
    E[k,j] = (Kcounter)/(Kcounter+1)*E[k,j]

#    F = -gam*ones(n,n)*sigmaK + (1-sigmaK)*F

end
end
end # end triple for loop


# Part 2: -E <= F and E <= F checking, the second set of constraints
for i = 1:n-1
  for j = i+1:n
    Kcounter += 1

    # Update E[j,i] and F[j,i] for the adjustments we haven't made since the last time we explicitly worked with them explicitly
    E[j,i] = (CE[j,i]+1)/Kcounter*E[j,i]
    CE[j,i] = Kcounter

    F[j,i] = -gam*(Kcounter-CF[j,i]-1)/(Kcounter) + (CF[j,i]+1)/Kcounter*F[j,i]
    CF[j,i] = Kcounter

    Eij = E[j,i]
    Fij = F[j,i]
    delta = -Eij - Fij
    if delta > epsi
      E[j,i] = Eij + delta/2
      F[j,i] = Fij + delta/2
    end

    # Now do the HLWB step with the steering sequence
    sigmaK = 1/(Kcounter +1)
    E[j,i] = (Kcounter)/(Kcounter+1)*E[j,i]
    F[j,i] = 1/(Kcounter+1)*(-gam) + (Kcounter/(Kcounter+1))*F[j,i]
#    F = -gam*ones(n,n)*sigmaK + (Kcounter)/(Kcounter+1)*F

end
end

for i = 1:n-1
  for j = i+1:n
    Kcounter += 1

    E[j,i] = (CE[j,i]+1)/Kcounter*E[j,i]
    CE[j,i] = Kcounter

    F[j,i] = -gam*(Kcounter-CF[j,i]-1)/(Kcounter) + (CF[j,i]+1)/Kcounter*F[j,i]
    CF[j,i] = Kcounter

    Eij = E[j,i]
    Fij = F[j,i]
    delta = Eij - Fij
    if delta > epsi
      E[j,i] = Eij - delta/2
      F[j,i] = Fij + delta/2
    end

    # Now do the HLWB step with the steering sequence

    sigmaK = 1/(Kcounter +1)
    E[j,i] = (1-sigmaK)*E[j,i]
    F[j,i] = 1/(Kcounter+1)*(-gam) + (Kcounter/(Kcounter+1))*F[j,i]
#    F = -gam*ones(n,n)*sigmaK + (1-sigmaK)*F

  end
end
# Done with the first time through the constraints

iter = 0
maxits = 100000
while true && iter < maxits
    tempE = E[:,:]

    # Start again with the triangle inequality constraints
    tic()
    # Part 1: Triangle inequality constraint checking
    @inbounds for i = 1:n-2
       for j = i+1:n-1
        for k = j+1:n

        Kcounter += 1
        sigmaK = 1/( Kcounter + 1)

        # Update Eij, Eik, Ejk for all the updates they missed of the form Eij = (1-sigmaK)*Eij
        E[j,i] = (CE[j,i]+1)/Kcounter*E[j,i]
        E[k,i] = (CE[k,i]+1)/Kcounter*E[k,i]
        E[k,j] = (CE[k,j]+1)/Kcounter*E[k,j]

        CE[j,i] = Kcounter
        CE[k,i] = Kcounter
        CE[k,j] = Kcounter

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


        # Now do the HLWB step with the steering sequence
        E[j,i] = (1-sigmaK)*E[j,i]
        E[k,i] = (1-sigmaK)*E[k,i]
        E[k,j] = (1-sigmaK)*E[k,j]

    #    F = -gam*ones(n,n)*sigmaK + (1-sigmaK)*F

    end
    end
    end # end triple for loop
    lasttime = toq()


    ### Now we check the non-triangle constraints
    # Part 2: -E <= F and E <= F checking, the second set of constraints
    for i = 1:n-1
      for j = i+1:n
        Kcounter += 1

        # Update E[j,i] and F[j,i] for the adjustments we haven't made since the last time we explicitly worked with them explicitly
        E[j,i] = (CE[j,i]+1)/Kcounter*E[j,i]
        CE[j,i] = Kcounter

        F[j,i] = -gam*(Kcounter-CF[j,i]-1)/(Kcounter) + (CF[j,i]+1)/Kcounter*F[j,i]
        CF[j,i] = Kcounter

        Eij = E[j,i]
        Fij = F[j,i]
        delta = -Eij - Fij
        if delta > epsi
          E[j,i] = Eij + delta/2
          F[j,i] = Fij + delta/2
        end

        # Now do the HLWB step with the steering sequence
        sigmaK = 1/(Kcounter +1)
        E[j,i] = (Kcounter)/(Kcounter+1)*E[j,i]
        F[j,i] = 1/(Kcounter+1)*(-gam) + (Kcounter/(Kcounter+1))*F[j,i]
    #    F = -gam*ones(n,n)*sigmaK + (Kcounter)/(Kcounter+1)*F

    end
    end

    for i = 1:n-1
      for j = i+1:n
        Kcounter += 1

        E[j,i] = (CE[j,i]+1)/Kcounter*E[j,i]
        CE[j,i] = Kcounter

        F[j,i] = -gam*(Kcounter-CF[j,i]-1)/(Kcounter) + (CF[j,i]+1)/Kcounter*F[j,i]
        CF[j,i] = Kcounter

        Eij = E[j,i]
        Fij = F[j,i]
        delta = Eij - Fij
        if delta > epsi
          E[j,i] = Eij - delta/2
          F[j,i] = Fij + delta/2
        end

        # Now do the HLWB step with the steering sequence

        sigmaK = 1/(Kcounter +1)
        E[j,i] = (1-sigmaK)*E[j,i]
        F[j,i] = 1/(Kcounter+1)*(-gam) + (Kcounter/(Kcounter+1))*F[j,i]
    #    F = -gam*ones(n,n)*sigmaK + (1-sigmaK)*F

      end
    end


    # Now before checking I must go through and make sure everything is updated
    # Kcounter is fixed, and I need to update any Eij such that CE[j,i] < Kcounter
    for i = 1:n-1
      for j = i+1:n
        if CE[j,i] < Kcounter
          E[j,i] = (CE[j,i]+1)/(Kcounter+1)*E[j,i]
          CE[j,i] = Kcounter
        end

        if CF[j,i] < Kcounter
          # update everything from start s = CF[j,i] to c = Kcounter
          s = CF[j,i]
          c = Kcounter
          F[j,i] = -gam*((c-s+1)/(c+1)) + s/(c+1)*F[j,i]
          CF[j,i] = Kcounter
        end
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
      tr = FullTriangleCheck(G)
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
