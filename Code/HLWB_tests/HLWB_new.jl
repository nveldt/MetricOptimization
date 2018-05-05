include("Tricon_Helper.jl")

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
function HLWB_new_TriangleFix(A::SparseMatrixCSC{Float64,Int64},D::Matrix{Int8},E::Matrix{Float64},F::Matrix{Float64},W::Matrix{Float64},Etol::Float64,TriTol::Float64,lam::Float64,filename::String,gam::Float64)
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

# Part 1: Triangle inequality constraint checking
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


# Part 2: -E <= F and E <= F checking, the second set of constraints
for i = 1:n-1
  for j = i+1:n

    Eij = E[j,i]
    Fij = F[j,i]
    delta = -Eij - Fij
    if delta > epsi
      E[j,i] = Eij + delta/2
      F[j,i] = Fij + delta/2
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
    end
  end
end
# Done with the first time through the constraints

# After going through and just projecting, now perform the update with steering sequence

# sigmaK = 1/(Kcounter+1)
# E = Kcounter/(Kcounter+1)*E
# F = -gam*ones(n,n)*sigmaK + (1-sigmaK)*F

con = 1/gam
sigmaK = con/(Kcounter+1)

for i = 1:n-1
  for j = i+1:n
    E[j,i] = (1-sigmaK)*E[j,i]
    F[j,i] = -gam*sigmaK + (1-sigmaK)*F[j,i]
  end
end


iter = 1
maxits = 1e5
while true && iter < maxits
    tempE = E[:,:]
    tempF = F[:,:]

    # Start again with the triangle inequality constraints
    tic()
    # Part 1: Triangle inequality constraint checking
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
    lasttime = toq()

    # Part 2: -E <= F and E <= F checking, the second set of constraints
    for i = 1:n-1
      for j = i+1:n

        Eij = E[j,i]
        Fij = F[j,i]
        delta = -Eij - Fij
        if delta > epsi
          E[j,i] = Eij + delta/2
          F[j,i] = Fij + delta/2
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
        end
      end
    end

    # Print and update of results for this loop
    iter += 1

    Kcounter = iter
    #con = 1/gam
    sigmaK = con/(Kcounter+1)

    for i = 1:n-1
      for j = i+1:n
        E[j,i] = (1-sigmaK)*E[j,i]
        F[j,i] = -gam*sigmaK + (1-sigmaK)*F[j,i]
      end
    end

    # E = Kcounter/(Kcounter+1)*E
    # F = -gam*ones(n,n)*sigmaK + (1-sigmaK)*F

    change = norm([vec(E-tempE); vec(F-tempF)])
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


return E+D, iter


end
