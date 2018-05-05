using MAT

include("Triangle_Fixing_LamCC2.jl")
include("Tri_Constraint_Library.jl")

mat = matread("/Users/nateveldt/GitHubRepos/TriConstraints/Code/KarateA.mat")
A = mat["A"]

mat = matread("/Users/nateveldt/GitHubRepos/TriConstraints/Code/graphs/SmallNets.mat")
Alesmis = mat["Alesmis"]
Afootball = mat["Afootball"]
Apolbooks = mat["Apolbooks"]
Adolphins = mat["Adolphins"]

A = Afootball

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
lam = .0001
W = lam*ones(n,n)
for i = 1:n-1
    for j = i+1:n
        if A[i,j] > .1
            W[j,i] = 1-lam
        end
    end
end

gam = 1000;
#gam = min(gam,n/2)

# Tolerance for how much we want matrix E to change before we declare convergence
Etol = .05

# Tolerance for how well the triangle constraints hold before declaring convergence
TriTol = .05

E = zeros(n,n);

# This is different. Instead of F = -gamma*W, we make it -gamma*ones,
# and then we use the W matrix when we are projecting
F = -gam*W;

# correction variables for the E <= F and -E <= F constraints
P = zeros(n,n)
Q = zeros(n,n)


tic()
M = TriangleFix(A,Dgraph,E,F,W,P,Q,Float64(gam),Float64(Etol),Float64(TriTol),lam,"TriFixOutput.txt",Float64(gam))
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
