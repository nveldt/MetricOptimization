using MAT

include("Serial_Han.jl")
include("Tri_Constraint_Library.jl")

mat = matread("/Users/nateveldt/GitHubRepos/TriConstraints/Code/KarateA.mat")
A = mat["A"]

mat = matread("/Users/nateveldt/GitHubRepos/TriConstraints/Code/graphs/SmallNets.mat")
Alesmis = mat["Alesmis"]
Afootball = mat["Afootball"]
Apolbooks = mat["Apolbooks"]
Adolphins = mat["Adolphins"]

A = Adolphins

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
lam = .3
W = lam*ones(n,n)
for i = 1:n-1
    for j = i+1:n
        if A[i,j] > .1
            W[j,i] = 1-lam
        end
    end
end

gam = 10;

# Tolerance for how much we want matrix E to change before we declare convergence
Etol = .1

# Tolerance for how well the triangle constraints hold before declaring convergence
TriTol = .1

E = zeros(n,n);

# This is different. Instead of F = -gamma*W, we make it -gamma*ones,
# and then we use the W matrix when we are projecting
F = -gam*ones(n,n);

# correction variables for the E <= F and -E <= F constraints
P = zeros(n,n)
Q = zeros(n,n)


tic()
M = TriangleFix(A,Dgraph,E,F,W,P,Q,Float64(gam),Float64(Etol),Float64(TriTol),lam,"TestHanOutput.txt",Float64(gam))
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

@show BoundFix
@show lccbound
