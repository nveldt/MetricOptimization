using MAT

include("Tri_Constraint_Library.jl")


name = "dolphinsA"
name = "lesmisA"
name = "BuckGraph"
name = "footballA"
name = "adjnounA"
name = "BuckGraph"
name = "KarateA"
name = "footballA"



mat = matread("graphs/"*name*".mat")
A = mat["A"]

gam = 75

n = 200
# mat = matread("graphs/A$n\_sparse.mat")
# A = mat["A"]

A = sparse(A)
n = size(A,1)

# tic()
# XLR, LP = LazyLeightonRao(A,true,gam)
# time1 = toc()
#
# @show LP

##
#
# xTx = 0.0
# for i = 1:n-1
# for j = i+1:n
#    xTx += XLR[j,i]^2
# end
# end
#
# TruePrimal = xTx/(2*gam)+ LR_obj(A,XLR)

##
# # tic()
# # Xfull, LPfull = LeightonRao(A)
# # timeFull = toc()
# #
# # @show timeFull, time1
# # @show LPfull, LP

##

include("LR_dualprimal2.jl")

TriTol = 1e-4
Xtol = 1e-1
gam = Float64(gam)
maxits = 50000
X, Finaltr, Finalobj, Finalits = Dykstra_LeightonRao_TFA(A,Xtol,TriTol,1.0,"boringfile",gam,maxits)

##
# include("LR_TFA1.jl")
# tol = .1
# Tritol = .00001
# lam = .1
#
# filename = "blah"
# gam = 100.0
# maxits = 10000
# tic()
# X, Finaltr, Finalobj, Finalits = Dykstra_LeightonRao_TFA(A,tol,Tritol,lam,filename,gam,maxits)
# tfatime = toq()
