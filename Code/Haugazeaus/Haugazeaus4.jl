include("Tricon_Helper.jl")

function cyclic_triangle_constraints!(D::Matrix{Int8},E::Matrix{Float64},W::Matrix{Float64})
n = size(D,1)
@inbounds for i = 1:n-2
    for j = i+1:n-1
        Wij = W[j,i]
        for k = j+1:n
            Wik = W[k,i]
            Wjk = W[k,j]
            denom = Wik*Wjk + Wij*Wik + Wij*Wjk
            # i,j,k here satisfies i < j < k
            # There are three ways to order these, each with a different constraint
            ABCproj!(D,E,i,j,Wij,Wik,Wjk,n,denom,k,i,k,j)
            ABCproj!(D,E,i,k,Wik,Wij,Wjk,n,denom,j,i,k,j)
            ABCproj!(D,E,j,k,Wjk,Wij,Wik,n,denom,j,i,k,i)
        end
    end
end
end

function ABCproj!(D::Matrix{Int8},E::Matrix{Float64},a::Int64,b::Int64,
  wab::Float64,wac::Float64,wbc::Float64,n::Int64,denom::Float64,maxAC::Int64,minAC::Int64,maxBC::Int64,minBC::Int64)

   # Get current iterates
   Eab = E[b,a]
   Eac = E[maxAC,minAC]
   Ebc = E[maxBC,minBC]
   Dab = D[b,a]
   Dac = D[maxAC,minAC]
   Dbc = D[maxBC,minBC]

   delta = Eab+Dab - Eac -Dac - Ebc - Dbc
   # We only change things if the constraint is violated
   if delta > 0
     E[b,a] -= (delta*wbc*wac/denom)
     E[maxAC,minAC] += (delta*wab*wbc/denom)
     E[maxBC,minBC] += (delta*wab*wac/denom)
  end
end


function doubleloop!(E::Matrix{Float64},F::Matrix{Float64},n::Int64)

  # Part 2: -E <= F and E <= F checking, the second set of constraints
  for i = 1:n-1
    for j = i+1:n

      Eij = E[j,i]
      Fij = F[j,i]
      delta = -Eij - Fij
      if delta > 0
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
      if delta > 0

        E[j,i] = Eij - delta/2
        F[j,i] = Fij + delta/2
    end
  end
end

end

function T1!(D::Matrix{Int8},E::Matrix{Float64},F::Matrix{Float64},W::Matrix{Float64})
    cyclic_triangle_constraints!(D,E,W)
end

function T2!(D::Matrix{Int8},E::Matrix{Float64},F::Matrix{Float64},W::Matrix{Float64})
    n = size(D,1)
    doubleloop!(E::Matrix{Float64},F::Matrix{Float64},n::Int64)
end

function Haugazeaus(A::SparseMatrixCSC{Float64,Int64},D::Matrix{Int8},E::Matrix{Float64},F::Matrix{Float64},W::Matrix{Float64},Etol::Float64,TriTol::Float64,lam::Float64,filename::String,gamma::Float64)
n = size(A,1)
@show filename
  open(filename, "w") do f
          write(f, "Writing output from calling the Triangle Fixing Algorithm\n")
          write(f, "Lambda = $lam, gamma = $gamma, Etol = $Etol, TriTol = $TriTol \n")
  end

# Tolerance for inequalty constraints
epsi = 1e-12
vareps = TriTol*3/4
tr = 0

ZK_E = 0

Elast = E[:,:]
Flast = F[:,:]

T1!(D,E,F,W)
HaugazeauStep!(E,Elast,F,Flast,gamma)

E_pm1_k = E[:,:]
F_pm1_k = F[:,:]

Elast = E[:,:]
Flast = F[:,:]

T2!(D,E,F,W)
HaugazeauStep!(E,Elast,F,Flast,gamma)

ZOE = E[:,:]
ZOF = F[:,:]
E_p_k = E[:,:]
F_p_k = F[:,:]
ze_0 = vec(E[:,:])
zf_0 = vec(F[:,:])
z0 = [ze_0;zf_0]

# Done with the first time through the constraints
iter = 0
maxits = 1e10
while true && iter < maxits


    tempE = E[:,:]
    tempF = F[:,:]

    Elast = E[:,:]
    Flast = F[:,:]
    # Now we map E to T(E) and similar for F
    tic()
    T1!(D,E,F,W)
    lasttime = toq()
    # Now perform the Haugazeus step
    HaugazeauStep!(E,Elast,F,Flast,gamma)

    # At this point E and F are the current iterate k+1
    Elast = E[:,:]
    Flast = F[:,:]

    E_pm1_kp1 = E[:,:]
    F_pm1_kp1 = F[:,:]

    # Ek and Fk are the current iterate
    T2!(D,E,F,W)
    # and now Ek and Fk have been updated to E and F

    # After another HaugazeauStep...
    HaugazeauStep!(E,Elast,F,Flast,gamma)
    # E and F are the new iterate

    E_p_kp1 = E[:,:]
    F_p_kp1 = F[:,:]

    # 3. Acceleration Part
    xpk_e = vec(E_p_k)
    xpk_f = vec(F_p_k)
    xpk = [xpk_e; xpk_f]

    xpkp1_e = vec(E_p_kp1)
    xpkp1_f = vec(F_p_kp1)
    xpkp1 = [xpkp1_e; xpkp1_f]

    xpm1kp1_e = vec(E_pm1_kp1)
    xpm1kp1_f = vec(F_pm1_kp1)
    xpm1kp1 = [xpm1kp1_e; xpm1kp1_f]

    vk_e = vec(E_p_k - E_pm1_k)
    vk_f = vec(F_p_k - F_pm1_k)
    vk = [vk_e; vk_f]

    vkp1_e = vec(E_pm1_kp1 - E_p_kp1)
    vkp1_f = vec(F_pm1_kp1 - F_p_kp1)
    vkp1 = [vkp1_e; vkp1_f]

    c1 = norm(vk)^2
    c2 = dot(vk,vkp1)

    if abs(c1+c2) > vareps
      ak = c1/(c1+c2)
      zk = xpk*(1-ak) + ak*(xpkp1)
      ZK_E = E_p_k + ak*(E_p_kp1 - E_p_k)
    else
      zk = .5*(xpkp1 + xpm1kp1)
      ZK_E = .5*(E_p_kp1 + E_pm1_kp1)
    end

    change = norm(zk - z0)

    z0 = zk

    # now k = k+1, so need to update the E and F vectors
    E_p_k = E_p_kp1[:,:]
    F_p_k = F_p_kp1[:,:]

    E_pm1_k = E_pm1_kp1[:,:]
    F_pm1_k = F_pm1_kp1[:,:]


    # Print and update of results for this loop
    iter += 1
    #change = norm([vec(E-tempE); vec(F-tempF)])
    tricheck = TriangleCheck(ZK_E+D,TriTol)
    objective = LPcc_obj(A,ZK_E+D,lam)
    ch = round(change,3)
    lt = round(lasttime,1)
    ob = round(objective,3)
    if iter%20 == 1
      G = ZK_E+D
      tr = FullTriangleCheck(G)
    end
    println("Iteration $iter, E changed by $ch, T-map: $lt, Obj: $ob, Last TriCheck = $tr")

    open(filename, "a") do f
      write(f,"Iteration $iter, E changed by $ch, T-map: $lt, Obj: $ob, Last TriCheck = $tr\n")
    end

    if change < Etol && tricheck
      break
    end

end #end while loop


return ZK_E+D

end

function HaugazeauStep!(E::Matrix{Float64},Ek::Matrix{Float64},F::Matrix{Float64},Fk::Matrix{Float64},gamma::Float64)
  Xe = 0
  Xf = 0
  ue = 0
  uf = 0
  vf = 0
  ve = 0
  for i = 1:n-1
    for j = i+1:n
      Xe += (-Ek[j,i])*(Ek[j,i] - E[j,i])
      Xf += (-gamma-Fk[j,i])*(Fk[j,i] - F[j,i])
      ue += Ek[j,i]^2
      uf += (gamma+F[j,i])^2
      ve += (Ek[j,i] - E[j,i])^2
      vf += (Fk[j,i] - F[j,i])^2
    end
  end
  v = vf+ve
  u = uf + ue
  X = Xe+Xf

  p = u*v - X*X

  if p == 0
    ## Case 1: Just take T(g) as the next iterate
    E = Ek[:,:]
    F = Fk[:,:]
    println("Something is weird if this evaluates to true")

  else # p > 0

    ## Case 2: Reset things
    if X*v >= p
      for i = 1:n-1
        for j = i+1:n
        E[j,i] = (1+X/v)*(E[j,i] - Ek[j,i])
        F[j,i] = -gamma + (1+X/v)*(F[j,i] - Fk[j,i])
        end
      end


    ## Case 3: Take a combination of new and old
    else
      # Main update

      for i = 1:n-1
        for j = i+1:n
          E[j,i] = (1-v/p*(X+u))*Ek[j,i] + u*v/p*E[j,i]
          F[j,i] = (1-v/p*(X+u))*Fk[j,i] + u*v/p*F[j,i] - v*X/p*gamma
        end
      end

    end
  end # end p == 0 conditional

end
