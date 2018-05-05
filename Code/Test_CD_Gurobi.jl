using MAT

#include("Tri_Constraint_Library.jl")
include("Cluster_Deletion_Library.jl")
name = "adjnounA"
name = "dolphinsA"
name = "lesmisA"
# name = "BuckGraph"
# name = "footballA"
# name = "KarateA"
name = "footballA"
mat = matread("graphs/"*name*".mat")
A = mat["A"]

n = 750
mat = matread("graphs/A$n\_sparse.mat")
A = mat["A"]

tic()
D, OneMinusD, cdBound = CD_relax_Gurobi(A)
timeCDrelax = toq()

tic()
Dold, cdBoundold = CD_relax_fullmemory(A,false)
timeCDrelaxOld = toq()

@show cdBound, cdBoundold, timeCDrelaxOld, timeCDrelax

cOld, scoreOld, approxOld = CDLP_dense_round(A,Dold,cdBoundold)
cOld2, scoreOld2, approxOld2 = CDLP_round(A,sparse(1-Dold),cdBoundold)
c2, score2, approx2 = CDLP_dense_round(A,full(1-OneMinusD),cdBound)
##
tic()
c, score, approx = CDLP_round(A,OneMinusD,cdBound)
first = toc()

tic()
c, score, approx = CDLP_round_slow(A,OneMinusD,cdBound)
second = toc()
