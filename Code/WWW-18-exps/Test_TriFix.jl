using MAT

include("TriFixSmallMemory.jl")
include("LamCC_LP.jl")

mat = matread("/Users/nateveldt/GitHubRepos/TriConstraints/Code/KarateA.mat")
A = mat["A"]


# Here I can check the solution for a three-node graph
# I'm trying to figure out whether the constraint x_{ij} <= 1 is needed
A = sparse(zeros(3,3))

A[1,2] = 1
A[2,3] = 1
A[2,3] = 1

A = A+A'

# sparse seems to NOT be the way to go
#mat = matread("../graphs/A500_sparse.mat")
#A = mat["A"]


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

lam = .5
W = lam*ones(n,n)

for i = 1:n-1
    for j = i+1:n
        if A[i,j] > .1
            W[j,i] = 1-lam
        end
    end
end
gam = max(10,1/(20*lam))
#gam = min(gam,n/2)

Etol = .1
TriTol = .1

E = zeros(n,n);
F = -gam*W;

# correction variables for the E <= F and -E <= F constraints
P = zeros(n,n)
Q = zeros(n,n)


tic()
M = TriangleFix(A,Dgraph,E,F,W,P,Q,Float64(gam),Float64(Etol),Float64(TriTol),lam,"TriFixOutput.txt",gam)
toc()

D = M+M'
#@show D[1:5,1:5]

# Check the bound
numedges = countnz(A)/2
BoundFix = sum((A[i,j]-lam)*D[i,j] for i=1:n-1 for j = i+1:n) + lam*(n*(n-1)/2 - numedges)

tic()
Dtrue, lccbound = FastLPlamCC(A,lam)
toc()

#@show Dtrue[1:5,1:5]

@show M+M'
@show Dtrue


# tic()
# Dtrue, lccbound = FastLPlamCCnobound(A,lam)
# toc()
