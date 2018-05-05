using MAT

include("New_Serial_Weighted_Tfix_safe.jl")
include("Tri_Constraint_Library.jl")
include("Experimental_HildrethsTFA.jl")

#include("Hildreths_sparsevec.jl")

#include("OldCode/Tfix_no_tri_corrections.jl")
mat = matread("KarateA.mat")
A = mat["A"]

mat = matread("graphs/SmallNets.mat")
Alesmis = mat["Alesmis"]
Afootball = mat["Afootball"]
Apolbooks = mat["Apolbooks"]
Adolphins = mat["Adolphins"]

A = Adolphins


n = 1000
mat = matread("graphs/A$n\_sparse.mat")
A = mat["A"]

n = size(A,1)
Dgraph = zeros(Int8,n,n)
for i = 1:n-1
    for j = i+1:n
        if A[i,j] < .1
            Dgraph[j,i] = 1
            #Dgraph[i,j] = 1
        end
    end
end

# The weights matrix will be the same as before, but it will be used differently
lam = .5
W = lam*ones(n,n)
for i = 1:n-1
    for j = i+1:n
        if A[i,j] > .1
            W[j,i] = 1-lam
        end
    end
end

gam = 10

# Tolerance for how much we want matrix E to change before we declare convergence
Etol = .25

# Tolerance for how well the triangle constraints hold before declaring convergence
TriTol = .1

E = zeros(n,n);

# This is different. Instead of F = -gamma*W, we make it -gamma*ones,
# and then we use the W matrix when we are projecting
F = -gam*ones(n,n);

# correction variables for the E <= F and -E <= F constraints
P = zeros(n,n)
Q = zeros(n,n)

#
# tic()
# M = TriangleFix(A,Dgraph,E,F,W,P,Q,Float64(gam),Float64(Etol),Float64(TriTol),lam,"TriFixOutput.txt",Float64(gam))
# trime = toc()
#
# D = M+M'
# #@show D[1:5,1:5]
#
# # Check the bound
# numedges = countnz(A)/2
# BoundFix = sum((A[i,j]-lam)*D[i,j] for i=1:n-1 for j = i+1:n) + lam*(n*(n-1)/2 - numedges)


# include("New_Serial_Weighted_Tfix_Safe.jl")
# E = zeros(n,n)
# F = -gam*ones(n,n)
# P = zeros(n,n)
# Q = zeros(n,n)
# tic()
# M = TriangleFix(A,Dgraph,E,F,W,P,Q,Float64(gam),Float64(Etol),Float64(TriTol),lam,"TriFixOutput.txt",Float64(gam))
# trimeOld = toc()
# D = M+M'
# BoundFixOld = sum((A[i,j]-lam)*D[i,j] for i=1:n-1 for j = i+1:n) + lam*(n*(n-1)/2 - numedges)


tic()
#Dtrue, lccbound = FastLPlamCC(A,lam)
LPtime = toc()


X = zeros(n,n)
D = zeros(Int8,n,n)
for i = 1:n-1
    for j = i+1:n
        if A[i,j] < .1
            #X[j,i] = 1
            D[j,i] = 1
        end
    end
end
W = lam*ones(n,n)
for i = 1:n-1
    for j = i+1:n
        if A[i,j] > .1
            W[j,i] = 1-lam
        end
    end
end
Y = -gam*ones(n,n)

## Matrices for double-loop correction terms (purposely dense)
P = zeros(n,n)
Q = zeros(n,n)
triplet_corrections = Dict{Int64,Float64}()
tic()
X, tri = Hildreth_lamCC_TFA(A,Etol,TriTol,lam,"hildreths",Float64(gam),1000,1.0,W,X,Y,P,Q,D,triplet_corrections)
hill = toc()
##
D = X+X'
HildrethObj = sum((A[i,j]-lam)*D[i,j] for i=1:n-1 for j = i+1:n) + lam*(n*(n-1)/2 - numedges)

@show HildrethObj, hill
@show BoundFix, trime
# @show BoundFixOld, trimeOld
#@show lccbound, LPtime
