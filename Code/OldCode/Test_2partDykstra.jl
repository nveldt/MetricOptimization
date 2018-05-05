using MAT

include("Two_part_Dykstra.jl")
include("Tri_Constraint_Library.jl")
mat = matread("/Users/nateveldt/GitHubRepos/Tri-Con-LP/Code/KarateA.mat")
A = mat["A"]


mat = matread("/Users/nateveldt/GitHubRepos/Tri-Con-LP/Code/graphs/SmallNets.mat")
Alesmis = mat["Alesmis"]
Afootball = mat["Afootball"]
Apolbooks = mat["Apolbooks"]
Adolphins = mat["Adolphins"]
# A = spones(Alesmis)
A = Adolphins

n = 500
# mat = matread("graphs/A$n\_sparse.mat")
# A = mat["A"]

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
Etol = .01
TriTol = .01
gam = 5

W = lam*ones(n,n)
for i = 1:n-1
    for j = i+1:n
        if A[i,j] > .1
            W[j,i] = 1-lam
        end
    end
end

# Initial "vector" X = [E;F].
E = zeros(n,n);
F = -gam*ones(n,n);

tic()
M = Dykstra_TriangleFix(A,Dgraph,E,F,W,Float64(Etol),Float64(TriTol),lam,"Dykstra_2constraints_Test.txt",Float64(gam))
Dykstra_time = toc()

D = M+M'

# Check the bound
numedges = countnz(A)/2
Dykstra_obj = sum((A[i,j]-lam)*D[i,j] for i=1:n-1 for j = i+1:n) + lam*(n*(n-1)/2 - numedges)

tic()
Dtrue, lccbound = FastLPlamCC(A,lam)
lazytime = toq()

#@show Dtrue[1:5,1:5]

@show lccbound, Dykstra_obj
@show lazytime, Dykstra_time
