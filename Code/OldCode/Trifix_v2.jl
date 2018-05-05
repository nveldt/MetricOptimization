# Code for Triangle Fixing Algorithms
#
# One based on Dykstra's method
# The other based on Bauchke's scheme

include("Tricon_Helper.jl")
# TriangleFix
#
# This solves the LambdaCC LP relaxation using the Triangle Fixing algorithm
# of Dhillon et al. for the L1 metric nearness problem.
#
# Basically it's a specific implementation of Dykstra's projection method
#
# D = 0,1 matrix indicating node distances--not satisfying triangle inequalities
# W = weights for each edge in the graph (based on what lambda is)
# gamma = parameter for turning the LP into a 2-norm problem (if this is not
#         set to a high enough value, this problem isn't actually equivalent to
#         the LP relaxation of correlation clustering.)
#
function Dykstra_TF(A::SparseMatrixCSC{Float64,Int64},D::Matrix{Int8},E::Matrix{Float64},F::Matrix{Float64},W::Matrix{Float64},P::Matrix{Float64},Q::Matrix{Float64},Etols::Array{Float64},TriTol::Float64,lam::Float64,file2::String,filename::String,gam::Float64,maxits::Int64)
n = size(A,1)

FinalTol = Etols[end]
open(filename, "w") do f
        write(f, "Output from Experiment on graph of size $n\n")
        write(f, "Lambda = $lam, gamma = $gam, Final Tolerance = $FinalTol, Triangle Constraint Tolerance = $TriTol \n")
end
open(file2, "w") do f
        write(f, "Intermediate convergence results\n Dykstra on graph of size $n\n")
        write(f, "Lambda = $lam, gamma = $gam, Final Tolerance = $FinalTol, Triangle Constraint Tolerance = $TriTol \n")
end

tr = 0

# Current Intermediate Tolerance
TolInd = 1
InTol = Etols[TolInd]
PrintNext = true

# We will fill this with triplets that need corrections in the first round
current_corrections = Vector{Tuple{Int64,Int64,Int64,Float64}}()

# Part 1: Triangle inequality constraint checking
FirstTripleLoop!(D,E,W,current_corrections)

# Part 2: -E <= F and E <= F checking, the second set of constraints
FirstDoubleLoop!(P,Q,E,F)

# In future rounds, we make corrections with now_corrections, while
# we prepare for new corrections with next_corrections
next_corrections = Vector{Tuple{Int64,Int64,Int64,Float64}}()

iter = 0
residual = 1
while true
    iter += 1
    tempE = E[:,:]
    tempF = F[:,:]

    # Start again with the triangle inequality constraints
    tic()
    TripleLoop!(D,E,W,current_corrections,next_corrections)
    lasttime = toq()

    # Get rid of old corrections, prepare space for new
    current_corrections = next_corrections
    next_corrections = Vector{Tuple{Int64,Int64,Int64,Float64}}()

    # Part 2: -E <= F and E <= F checking
    DoubleLoop!(P,Q,E,F)

    # Print and update of results for this loop
    residual = norm([vec(E-tempE); vec(F-tempF)])
    tricheck = TriangleCheck(E+D,TriTol)
    objective = LPcc_obj(A,E+D,lam)
    ch = round(residual,3)
    lt = round(lasttime,1)
    ob = round(objective,3)
    if iter%20 == 1
      tr = FullTriangleCheck(E+D)
    end

    println("Dykstra: Iter $iter, Residual = $ch, 3Loop Time = $lt, Obj: $ob, Last TriCheck = $tr")

    open(filename, "a") do f
      write(f,"Iter $iter, Residual = $ch, 3Loop Time = $lt, Obj: $ob, Last TriCheck = $tr\n")
    end

    # Check for final convergence
    if residual < FinalTol && tricheck
      open(file2, "a") do f
              write(f, "Final: Iter $iter, Residual = $ch < Tol = $FinalTol, LP objective = $ob, Triangle Constraint Satisfaction: $tr\n")
      end
      break
    end

    # Check for "convergence" to within an intermediate value
    if residual < InTol && PrintNext && TriangleCheck(E+D,.5)
      tr = FullTriangleCheck(E+D)
      open(file2, "a") do f
              write(f, "Iter $iter, Residual = $ch < Tol  = $InTol,  LP objective = $ob, Triangle Constraint Satisfaction: $tr\n")
      end
      PrintNext = false
      if TolInd < length(Etols)
         TolInd += 1
         InTol = Etols[TolInd]
         PrintNext = true
      end
    end

    # Check whether we reached maxits
    if iter >= maxits
        open(filename, "a") do f
              write(f,"Maximum number of iterations reached. Terminating algorithm.\n")
        end
        tr = FullTriangleCheck(E+D)
        open(file2, "a") do f
                write(f, "Terminated after $maxits iterations. Final residual = $ch: LP objective = $ob, Triangle Constraint Satisfaction: $tr\n")
        end
        break
    end

end #end while loop

# Return the objective value, the triangle constraint satisfaction, the number of iterations, and the last residual
Finaltr = FullTriangleCheck(E+D)
Finalobj = LPcc_obj(A,E+D,lam)
Finalits = iter
return Finalobj, Finalits, Finaltr, residual

end

# Now the Bauschke Projection version Triangle Fixing Algorithm
function Bauschke_TF(A::SparseMatrixCSC{Float64,Int64},D::Matrix{Int8},E::Matrix{Float64},F::Matrix{Float64},W::Matrix{Float64},Etols::Array{Float64},TriTol::Float64,lam::Float64,file2::String,filename::String,gam::Float64,sig::Float64,vareps::Float64,maxits::Int64)
n = size(A,1)
FinalTol = Etols[end]
open(filename, "w") do f
        write(f, "Output from Experiment on graph of size $n\n")
        write(f, "Lambda = $lam, gamma = $gam, Final Tolerance = $FinalTol, Triangle Constraint Tolerance = $TriTol \n")
end
open(file2, "w") do f
        write(f, "Intermediate convergence results\n Bauschke on Graph of size $n\n")
        write(f, "Lambda = $lam, gamma = $gam, Final Tolerance = $FinalTol, Triangle Constraint Tolerance = $TriTol \n")
end


tr = 1

# Current Intermediate Tolerance
TolInd = 1
InTol = Etols[TolInd]
PrintNext = true

ZK_E = 0

# 1. Triangle constraints
TripleLoop_NoMemory!(D,E,W)

K = 1; sigmaK = sig/(1+K)

Bauschke_Update!(E,F,sigmaK,n,gam)

E_pm1_k = E[:,:]
F_pm1_k = F[:,:]

# 2. Non-triangle constraints
DoubleLoop_NoMemory!(E,F)

K = 2
sigmaK = sig/(1+K)
Bauschke_Update!(E,F,sigmaK,n,gam)

ZOE = E[:,:]
ZOF = F[:,:]
E_p_k = E[:,:]
F_p_k = F[:,:]
ze_0 = vec(E[:,:])
zf_0 = vec(F[:,:])
z0 = [ze_0;zf_0]
residual = 1
iter = 1
while true
    iter += 1
    tempE = E[:,:]
    tempF = F[:,:]

    # 1. Triangle Constraints
    tic()
    TripleLoop_NoMemory!(D,E,W)
    lasttime = toq()

    K += 1
    sigmaK = sig/(1+K)
    Bauschke_Update!(E,F,sigmaK,n,gam)

    E_pm1_kp1 = E[:,:]
    F_pm1_kp1 = F[:,:]

    # 2. Non-Triangle Constraints
    DoubleLoop_NoMemory!(E,F)

    K += 1
    sigmaK = sig/(1+K)
    Bauschke_Update!(E,F,sigmaK,n,gam)

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

    residual = norm(zk - z0)

    z0 = zk

    # now k = k+1, so need to update the E and F vectors
    E_p_k = E_p_kp1[:,:]
    F_p_k = F_p_kp1[:,:]

    E_pm1_k = E_pm1_kp1[:,:]
    F_pm1_k = F_pm1_kp1[:,:]

    #change = norm([vec(E-tempE); vec(F-tempF)])
    tricheck = TriangleCheck(ZK_E+D,TriTol)
    objective = LPcc_obj(A,ZK_E+D,lam)
    ch = round(residual,3)
    lt = round(lasttime,1)
    ob = round(objective,3)
    if iter%20 == 1
        tr = FullTriangleCheck(ZK_E+D)
    end
    println("Bauschke: Iter $iter, Residual = $ch, 3Loop: $lt, Obj: $ob, Last TriCheck = $tr")

    open(filename, "a") do f
      write(f,"Bauschke: Iter $iter, Residual = $ch, 3Loop: $lt, Obj: $ob, Last TriCheck = $tr\n")
    end

    # Check for final convergence
    if residual < FinalTol && tricheck
      open(file2, "a") do f
              write(f, "Final: Iter $iter, Residual = $ch < Tol = $FinalTol, LP objective = $ob, Triangle Constraint Satisfaction: $tr\n")
      end
      break
    end

    # Check for "convergence" to within an intermediate value
    if residual < InTol && PrintNext && TriangleCheck(ZK_E+D,.5)
      tr = FullTriangleCheck(E+D)
      open(file2, "a") do f
              write(f, "Iter $iter, Residual = $ch < Tol  = $InTol,  LP objective = $ob, Triangle Constraint Satisfaction: $tr\n")
      end
      PrintNext = false
      if TolInd < length(Etols)
         TolInd += 1
         InTol = Etols[TolInd]
         PrintNext = true
      end
    end

    # Check whether we reached maxits
    if iter >= maxits
        open(filename, "a") do f
              write(f,"Maximum number of iterations reached. Terminating algorithm.\n")
        end
        tr = FullTriangleCheck(ZK_E+D)
        open(file2, "a") do f
                write(f, "Terminated after $maxits iterations. Final residual = $ch: LP objective = $ob, Triangle Constraint Satisfaction: $tr\n")
        end
        break
    end

 end #end while loop

 # Return the objective value, the triangle constraint satisfaction, the number of iterations, and the last residual
 Finaltr = FullTriangleCheck(ZK_E+D)
 Finalobj = LPcc_obj(A,ZK_E+D,lam)
 Finalits = iter
 return Finalobj, Finalits, Finaltr, residual

 end
