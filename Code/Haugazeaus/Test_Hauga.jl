using MAT

include("Haugazeaus2.jl")
include("New_Serial_Weighted_Tfix.jl")
include("Tri_Constraint_Library.jl")

mat = matread("KarateA.mat")
A = mat["A"]


mat = matread("graphs/SmallNets.mat")
Alesmis = mat["Alesmis"]
Afootball = mat["Afootball"]
Apolbooks = mat["Apolbooks"]
Adolphins = mat["Adolphins"]
#
# A = spones(Alesmis)
#A = Adolphins

n = 300
mat = matread("graphs/A$n\_sparse.mat")
A = mat["A"]

n = size(A,1)
Dgraph = zeros(Int8,n,n)
for i = 1:n-1
    for j = i+1:n
        if A[i,j] < .1
            Dgraph[j,i] = 1
        end
    end
end

# The weights matrix will be the same as before, but it will be used differently

lam = .001
tic()
# Dtrue, lccbound = FastLPlamCC(A,lam)
toc()

W = lam*ones(n,n)
for i = 1:n-1
    for j = i+1:n
        if A[i,j] > .1
            W[j,i] = 1-lam
        end
    end
end
gam = 5
#gam = min(gam,n/2)

# Tolerance for how much we want matrix E to change before we declare convergence
Etol = .1

# Tolerance for how well the triangle constraints hold before declaring convergence
TriTol = .1

## Old way


E = zeros(n,n);
F = -gam*ones(n,n);

# correction variables for the E <= F and -E <= F constraints
P = zeros(n,n)
Q = zeros(n,n)

tic()
M = TriangleFix(A,Dgraph,E,F,W,P,Q,Float64(gam),Float64(Etol),Float64(TriTol),lam,"TriFixOutput.txt",Float64(gam))
Dtime = toc()

D = M+M'
numedges = countnz(A)/2
DykstraObj = sum((A[i,j]-lam)*D[i,j] for i=1:n-1 for j = i+1:n) + lam*(n*(n-1)/2 - numedges)
#
# include("Haugazeaus4.jl")
# ## Haugazeaus 3
# E = zeros(n,n)
# F = -gam*tril(ones(n,n))
# tic()
# M = Haugazeaus(A,Dgraph,E,F,W,Float64(Etol),Float64(TriTol),lam,"Haugazeaus.txt",Float64(gam))
# Haug3 = toc()
# D = M+M'
# HaugObj3 = sum((A[i,j]-lam)*D[i,j] for i=1:n-1 for j = i+1:n) + lam*(n*(n-1)/2 - numedges)

#
#
# ## Haugazeaus 2
# E = zeros(n,n)
# F = -gam*tril(ones(n,n))
# tic()
# M = Haugazeaus(A,Dgraph,E,F,W,Float64(Etol),Float64(TriTol),lam,"Haugazeaus.txt",Float64(gam))
# Haug2 = toc()
# D = M+M'
# HaugObj2 = sum((A[i,j]-lam)*D[i,j] for i=1:n-1 for j = i+1:n) + lam*(n*(n-1)/2 - numedges)
#



#@show lccbound
# @show HaugObj3, Haug3
@show DykstraObj, Dtime
