# This is code specifically for running the Triangle Fixing Algorithm
# for finding LP lower bounds for LambdaCC on the ca-GrQc graph.
# (See experiment 2 in Veldt et al. WWW 2018)

using MAT

function run_cagrqc(filename::String,lam::Float64,Etol::Float64,TriTol::Float64,gam::Float64)

include("../Weighted_Serial_Tfix.jl")

# path = "/homes/lveldt/GitHubRepos/"
file = "Acagrqc_conn"
mat = matread(file*".mat")

A = mat["A"]
n = size(A,1)

# Set up initial matrix that we want to approximate with a matrix
# satisfying metric constraints
Dgraph = zeros(Int8,n,n)
for i = 1:n-1
    for j = i+1:n
        if A[i,j] < .1
            Dgraph[j,i] = 1
            #Dgraph[i,j] = 1
        end
    end
end

# Set up weights matrix for the instance of LambdaCC
W = lam*ones(n,n)

for i = 1:n-1
    for j = i+1:n
        if A[i,j] > .1
            W[j,i] = 1-lam
        end
    end
end

# Set up variables for the algorithm
E = zeros(n,n);
F = -gam*ones(n,n);

# correction variables for the E <= F and -E <= F constraints
P = zeros(n,n)
Q = zeros(n,n)

tic()
M = TriangleFix(A,Dgraph,E,F,W,P,Q,gam,Etol,TriTol,lam,filename,gam)
toc()
@show norm(vec(M))

D = M+M'
# Check the bound
numedges = countnz(A)/2
BoundFix = sum((A[i,j]-lam)*D[i,j] for i=1:n-1 for j = i+1:n) + lam*(n*(n-1)/2 - numedges)

@show BoundFix

end
