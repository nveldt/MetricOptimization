using MAT

include("../Algs/Trifix_v2.jl")
include("Tri_Constraint_Library.jl")
include("Hildreth_Xvar/Experimental_HildrethsTFA.jl")


name = "adjnounA"
name = "dolphinsA"
name = "lesmisA"
name = "BuckGraph"
name = "footballA"
name = "KarateA"
#name = "footballA"
mat = matread("graphs/"*name*".mat")
A = mat["A"]

n = 200
mat = matread("graphs/A$n\_sparse.mat")
A = mat["A"]
#
# n = 200
# A = zeros(n,n)
#
# for i = 1:n
#     for j = i+1:n
#         if rand(1)[1] > .5
#             A[i,j] = 1
#             A[j,i] = 1
#         end
#     end
# end
# A = sparse(A)

n = size(A,1)
D = zeros(Int8,n,n)
for i = 1:n-1
    for j = i+1:n
        if A[i,j] < .1
            D[j,i] = 1
        end
    end
end

# The weights matrix will be the same as before, but it will be used differently
lam = .1
W = lam*ones(n,n)
for i = 1:n-1
    for j = i+1:n
        if A[i,j] > .1
            W[j,i] = 1-lam
        end
    end
end

gam = 100

# Objective function tolerance and constraint tolerance
tol = .01
ConTol = .01

## Older, working code for Dykstra's version of the TFA
# E = zeros(n,n);
# F = -gam*ones(n,n);
# P = zeros(n,n)
# Q = zeros(n,n)
# #
# tic()
# Finalobj, Finalits, Finaltr, residual = Dykstra_TF(A,D,E,F,W,P,Q,[tol, tol/2],ConTol,lam,"short_file","Dykstra_output",Float64(gam),1000)
# Basetime = toq()

## Gurobi, lazy constraints method
# tic()
# Dtrue, lccbound = FastLPlamCC(A,lam)
# LPtime = toc()

##
include("Conv-Criteria/Dykstra_CCsafe.jl")
# Approximation factor tolerance
ApproxTol = 1.0000000001
tol = 1e-10
tic()
Xcc, cctr, ccobj, ccits, Ratio = Dykstra_lamCC_TFA(A,ApproxTol,tol,lam,"DykstraLamCCoutput",Float64(gam),10000)
newtime = toq()

# @show lccbound
# @show Finalobj, Finalits, Finaltr, residual
# @show cctr2, ccobj2, ccits2, Ratio2
# @show cctr, ccobj, ccits, Ratio
