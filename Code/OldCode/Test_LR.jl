include("Tri_Constraint_Library.jl")

using MAT


name = "dolphinsA"
name = "lesmisA"
name = "BuckGraph"
name = "footballA"
name = "adjnounA"
name = "BuckGraph"
name = "KarateA"
name = "footballA"

name = "jazzA"
#name = "BuckGraph"
#name = "Netscience"
mat = matread("graphs/"*name*".mat")
A = mat["A"]
A = sparse(mat["A"])

# n = 800
# mat = matread("graphs/A$n\_sparse.mat")
# A = mat["A"]

n = size(A,1)
using MatrixNetworks

##
# Setting Gamma may be helped by computing a different relaxation
# Actually, this is the second smallest eigval of the normalized laplacian,
# but the theory I developed is for the standard Laplacian L = D - A
#(v,lam2) = fiedler_vector(A)

d = sum(A,2)
D = zeros(n,n)
for i = 1:n
    D[i,i] = d[i]
end
L = D-A
E = eigvals(L)
lam2 = E[2]
gamma = n/(2*lam2)
##
# A = sparse(A)
# tic()
# XLR, LP = LazyLeightonRao(A,false,75)
# time1 = toc()
#
# @show LP

#include("LR_dualprimal2.jl")

# The dictionary version seems to work, at least the dual is increasing
include("CodePackage/LR_Dict.jl")
TriTol = 1e-5
GapTol = 1e-4
gam = Float64(1)
maxits = 5000
lami = lam2/n*2
#lami = 1.0
# X, Finaltr, Finalobj, Finalits = Dykstra_LeightonRao_TFA(A,GapTol,TriTol,lami,"LR1",gam,maxits)

##

include("CodePackage/TFA_LeightonRao.jl")

tic()
X, Finaltr, Finalobj, Finalits = Dykstra_LeightonRao_TFA(A,GapTol,TriTol,lami,"LR2",gam,maxits,10,1e-19)
DykstraTime = toq()

@show Finalobj, DykstraTime

##
# @show norm(vec(triu(X)))
# @show LP
# #
# tic()
# Xfull, LPfull = LeightonRao(A)
# timeFull = toc()
# #
# # @show timeFull, time1
# # @show LPfull, LP
##

##
# ## Let's try the Triangle Fixing Leighton Rao
# include("LR_Bauschke.jl")
# #include("LR_TFA1.jl")
# tol = .1
# Tritol = .00001
# lam = 1.0
#
# filename = "BauschkeLeightonRaoOutput"
# gam = 100.0
# maxits = 10000
# tic()
# X, Finaltr, Finalobj, Finalits = Dykstra_LeightonRao_TFA(A,tol,Tritol,lam,filename,gam,maxits)
# tfatime = toq()
#
# @show tfatime, Finalobj
