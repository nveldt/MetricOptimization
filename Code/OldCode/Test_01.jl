using MAT

include("Trifix_Library.jl")
include("Tri_Constraint_Library.jl")
mat = matread("KarateA.mat")
A = mat["A"]

mat = matread("WWW-18-exps/SmallNets.mat")
Alesmis = mat["Alesmis"]
Afootball = mat["Afootball"]
Apolbooks = mat["Apolbooks"]
Adolphins = mat["Adolphins"]
#
# A = spones(Alesmis)
A = Adolphins

n = 500
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
Etol = .01
TriTol = .05

E = zeros(n,n)
F = -gam*ones(n,n)
tic()
Mold = Dykstra_TriangleFix01bounds(A,Dgraph,E,F,W,Float64(Etol),Float64(TriTol),lam,"WeightTest.txt",Float64(gam))
Dykstra01 = toc()

E = zeros(n,n)
F = -gam*ones(n,n)
tic()
Mold = Dykstra_TriangleFix01bounds(A,Dgraph,E,F,W,Float64(Etol),Float64(TriTol),lam,"WeightTest.txt",Float64(gam))
Dykstra01 = toc()

D = Mold+Mold'

numedges = countnz(A)/2
ZeroOne = sum((A[i,j]-lam)*D[i,j] for i=1:n-1 for j = i+1:n) + lam*(n*(n-1)/2 - numedges)

## Run Old version twice
E = zeros(n,n);
F = -gam*ones(n,n);
P = zeros(n,n)
Q = zeros(n,n)
tic()
M = Dykstra_TriangleFix(A,Dgraph,E,F,W,P,Q,Float64(Etol),Float64(TriTol),lam,"WeightTest.txt",Float64(gam))
Dykstra = toc()
D = M+M'

E = zeros(n,n);
F = -gam*ones(n,n);
P = zeros(n,n)
Q = zeros(n,n)
tic()
Mold = Dykstra_TriangleFix(A,Dgraph,E,F,W,P,Q,Float64(Etol),Float64(TriTol),lam,"WeightTest.txt",Float64(gam))
Dykstra = toc()

D = Mold+Mold'
#@show D[1:5,1:5]

# Check the bound
numedges = countnz(A)/2
New = sum((A[i,j]-lam)*D[i,j] for i=1:n-1 for j = i+1:n) + lam*(n*(n-1)/2 - numedges)


# tic()
# Dtrue, lccbound = FastLPlamCC(A,lam)
# toc()
#
# #@show Dtrue[1:5,1:5]
#
# @show lccbound
@show ZeroOne, New
@show Dykstra01, Dykstra
