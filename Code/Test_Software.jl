using MAT
using Plots

include("../Algs/Trifix_v2.jl")
include("Tri_Constraint_Library.jl")
include("CodePackage/TFA_Projections_Library.jl")
include("CodePackage/TFA_Library_Helper.jl")

name = "adjnounA"
name = "dolphinsA"
name = "lesmisA"
name = "BuckGraph"
name = "footballA"
#name = "KarateA"
#name = "footballA"
mat = matread("graphs/"*name*".mat")
A = mat["A"]

n = 500
mat = matread("graphs/A$n\_sparse.mat")
A = mat["A"]

##
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


# The weights matrix will be the same as before, but it will be used differently
lam = .5

n = size(A,1)

D = zeros(Float64,n,n)
W = (1-lam)*ones(n,n)
for i = 1:n-1
    for j = i+1:n
        if A[i,j] < .1
            D[j,i] = 1
            W[j,i] = lam
        end
    end
end

# Gurobi, lazy constraints method
# tic()
# Dtrue, lccbound = FastLPlamCC(A,lam)
# LPtime = toc()
#
# xWx = Wnorm(W,Dtrue-D,abs.(Dtrue-D))
#
# TrueDual = lccbound + xWx/(2*gam)

##
# Approximation factor tolerance
gam = 20
GapTol = 1e-4
ConTol = 1e-5
MaxIts = Int64(500)
statusFreq = 10
stagTol = 1e-9
tic()
X, LPobjs, primals, duals, gaps, ConViolation, FinalCon, Finalits, Finalobj, R = Dykstra_lamCC_TFA(A,
GapTol,ConTol,lam,"DykstraLamCCoutput",Float64(gam),MaxIts,statusFreq,stagTol)
newtime = toq()

approx = (1+1/gam)/(1+R)

# plot(primals)
# plot!(duals)
# plot!(LPobjs)
