include("Tricon_Helper.jl")

# A two-step approach to HLWB, where we alternate between iteratively projecting
# onto triangle constraints, and then alternate between projecting onto
# non-triangle constraints one at a time.
function HLWB_TwoStep(A::SparseMatrixCSC{Float64,Int64},D::Matrix{Int8},E::Matrix{Float64},F::Matrix{Float64},W::Matrix{Float64},Etol::Float64,TriTol::Float64,lam::Float64,filename::String,gam::Float64,sig::Float64)
n = size(A,1)
@show filename
  open(filename, "w") do f
          write(f, "Writing output from calling the Triangle Fixing Algorithm, HLWB naive style\n")
          write(f, "Lambda = $lam, gamma = $gam, Etol = $Etol, TriTol = $TriTol \n")
  end

# Tolerance for inequality constraints
tr =  1

TripleLoop_NoMemory!(D,E,W)

K = 1
sigmaK = sig/(1+K)

for i = 1:n-1
  for j = i+1:n
    E[j,i] = (1-sigmaK)*E[j,i]
    F[j,i] = -gam*sigmaK + (1-sigmaK)*F[j,i]
  end
end

DoubleLoop_NoMemory!(E,F)

K = 2
sigmaK = sig/(1+K)
for i = 1:n-1
  for j = i+1:n
    E[j,i] = (1-sigmaK)*E[j,i]
    F[j,i] = -gam*sigmaK + (1-sigmaK)*F[j,i]
  end
end

iter = 1
maxits = 1e5
while true && iter < maxits
    iter += 1
    tempE = E[:,:]
    tempF = F[:,:]

    # 1. Triangle Constraints
    tic()
    TripleLoop_NoMemory!(D,E,W)
    lasttime = toq()

    K += 1
    sigmaK = sig/(1+K)


    for i = 1:n-1
      for j = i+1:n
        E[j,i] = (1-sigmaK)*E[j,i]
        F[j,i] = -gam*sigmaK + (1-sigmaK)*F[j,i]
      end
    end

    # 2. Non-Triangle Constraints
    DoubleLoop_NoMemory!(E,F)

    K += 1
    sigmaK = sig/(1+K)
    for i = 1:n-1
      for j = i+1:n
        E[j,i] = (1-sigmaK)*E[j,i]
        F[j,i] = -gam*sigmaK + (1-sigmaK)*F[j,i]
      end
    end

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


function HLWB_OneStep(A::SparseMatrixCSC{Float64,Int64},D::Matrix{Int8},E::Matrix{Float64},F::Matrix{Float64},W::Matrix{Float64},Etol::Float64,TriTol::Float64,lam::Float64,filename::String,gam::Float64,sig::Float64)
  n = size(A,1)
@show filename
  open(filename, "w") do f
          write(f, "Writing output from calling the Triangle Fixing Algorithm, HLWB naive style\n")
          write(f, "Lambda = $lam, gamma = $gam, Etol = $Etol, TriTol = $TriTol \n")
  end

# Tolerance for inequality constraints
tr =  1
epsi  = 0
TripleLoop_NoMemory!(D,E,W)
DoubleLoop_NoMemory!(E,F)

K = 1
sigmaK = sig/(1+K)

for i = 1:n-1
  for j = i+1:n
    E[j,i] = (1-sigmaK)*E[j,i]
    F[j,i] = -gam*sigmaK + (1-sigmaK)*F[j,i]
  end
end

iter = 1
maxits = 1e5
while true && iter < maxits
    iter += 1
    tempE = E[:,:]
    tempF = F[:,:]

    # 1. Triangle Constraints
    tic()
    TripleLoop_NoMemory!(D,E,W)
    lasttime = toq()

    DoubleLoop_NoMemory!(E,F)

    K += 1
    sigmaK = sig/(1+K)
    for i = 1:n-1
      for j = i+1:n
        E[j,i] = (1-sigmaK)*E[j,i]
        F[j,i] = -gam*sigmaK + (1-sigmaK)*F[j,i]
      end
    end

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


function HLWB_TwoStep_Accelerated(A::SparseMatrixCSC{Float64,Int64},D::Matrix{Int8},E::Matrix{Float64},F::Matrix{Float64},W::Matrix{Float64},Etol::Float64,TriTol::Float64,lam::Float64,filename::String,gam::Float64,sig::Float64,vareps::Float64)
n = size(A,1)
@show filename
  open(filename, "w") do f
          write(f, "Writing output from calling the Triangle Fixing Algorithm, HLWB naive style\n")
          write(f, "Lambda = $lam, gamma = $gam, Etol = $Etol, TriTol = $TriTol \n")
  end
tr =  1

ZK_E = 0

# 1. Triangle constraints
TripleLoop_NoMemory!(D,E,W)

K = 1
sigmaK = sig/(1+K)

for i = 1:n-1
  for j = i+1:n
    E[j,i] = (1-sigmaK)*E[j,i]
    F[j,i] = -gam*sigmaK + (1-sigmaK)*F[j,i]
  end
end

E_pm1_k = E[:,:]
F_pm1_k = F[:,:]

# 2. Non-triangle constraints
DoubleLoop_NoMemory!(E,F)

K = 2
sigmaK = sig/(1+K)
for i = 1:n-1
  for j = i+1:n
    E[j,i] = (1-sigmaK)*E[j,i]
    F[j,i] = -gam*sigmaK + (1-sigmaK)*F[j,i]
  end
end

ZOE = E[:,:]
ZOF = F[:,:]
E_p_k = E[:,:]
F_p_k = F[:,:]
ze_0 = vec(E[:,:])
zf_0 = vec(F[:,:])
z0 = [ze_0;zf_0]

iter = 1
maxits = 1e5
while true && iter < maxits
    iter += 1
    tempE = E[:,:]
    tempF = F[:,:]

    # 1. Triangle Constraints
    tic()
    TripleLoop_NoMemory!(D,E,W)
    lasttime = toq()

    K += 1
    sigmaK = sig/(1+K)

    for i = 1:n-1
      for j = i+1:n
        E[j,i] = (1-sigmaK)*E[j,i]
        F[j,i] = -gam*sigmaK + (1-sigmaK)*F[j,i]
      end
    end

    E_pm1_kp1 = E[:,:]
    F_pm1_kp1 = F[:,:]

    # 2. Non-Triangle Constraints
    DoubleLoop_NoMemory!(E,F)

    K += 1
    sigmaK = sig/(1+K)
    for i = 1:n-1
      for j = i+1:n
        E[j,i] = (1-sigmaK)*E[j,i]
        F[j,i] = -gam*sigmaK + (1-sigmaK)*F[j,i]
      end
    end

    # Update Current x_k
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
    println("Iteration $iter, E changed by $ch, 3Loop: $lt, Obj: $ob, Last TriCheck = $tr")

    open(filename, "a") do f
      write(f,"Iteration $iter, E changed by $ch, 3Loop: $lt, Obj: $ob, Last TriCheck = $tr\n")
    end

    if change < Etol && tricheck
      break
    end

end #end while loop


return ZK_E+D, iter


end
