using MAT

include("Trifix_Library.jl")
include("Tri_Constraint_Library.jl")
include("HLWB_new.jl")


# Super easy test case
mat = matread("KarateA.mat")
A = mat["A"]

# Slightly Larger test cases (be sure to make unweighted)
mat = matread("/Users/nateveldt/GitHubRepos/Tri-Con-LP/Code/graphs/SmallNets.mat")
Alesmis = mat["Alesmis"]
Afootball = mat["Afootball"]
Apolbooks = mat["Apolbooks"]
Adolphins = mat["Adolphins"]

# A = spones(Alesmis)
A = Adolphins

# Larger Test Cases
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

# Weights matrix W
lam = .5
gam = 50

Etol = .1
TriTol = .1

W = lam*ones(n,n)
for i = 1:n-1
    for j = i+1:n
        if A[i,j] > .1
            W[j,i] = 1-lam
        end
    end
end

E = zeros(n,n);
F = -gam*ones(n,n);

# correction variables for the E <= F and -E <= F constraints
P = zeros(n,n)
Q = zeros(n,n)

tic()
Dyout = Dykstra_TriangleFix(A,Dgraph,E,F,W,P,Q,Float64(Etol),Float64(TriTol),lam,"Dykstra_Tfix.txt",Float64(gam))
Dykstra_time = toc()

DyIt = Dyout[2]
M = Dyout[1]

D = M+M'

# Check the bound
numedges = countnz(A)/2
Dykstra_obj = sum((A[i,j]-lam)*D[i,j] for i=1:n-1 for j = i+1:n) + lam*(n*(n-1)/2 - numedges)

include("Twostep_hwlb.jl")

# 2. The One-Step HLWB
E = zeros(n,n)
F = -gam*ones(n,n)

Etol = .01
TriTol = .01
tic()
Hout =  HLWB_OneStep(A,Dgraph,E,F,W,Float64(Etol),Float64(TriTol),lam,"HLWB_new_Tfix.txt",Float64(gam),Float64(1/gam))
HLWB_new_time = toc()
hIt = Hout[2]
M = Hout[1]
D2 = M+M'
HLWB_new_obj = sum((A[i,j]-lam)*D2[i,j] for i=1:n-1 for j = i+1:n) + lam*(n*(n-1)/2 - numedges)

# 3. The Two-Step HLWB
E = zeros(n,n)
F = -gam*ones(n,n)
tic(); Hout2 = HLWB_TwoStep(A,Dgraph,E,F,W,Float64(Etol),Float64(TriTol),lam,"HLWB_new_Tfix.txt",Float64(gam),Float64(1/gam)); HLWB_new_time2 = toc()
hIt2 = Hout2[2]
M = Hout2[1]
D2 = M+M'
HLWB_new_obj2= sum((A[i,j]-lam)*D2[i,j] for i=1:n-1 for j = i+1:n) + lam*(n*(n-1)/2 - numedges)

# 4. The Accelerated Two-Step HLWB
E = zeros(n,n)
F = -gam*ones(n,n)
vareps = Float64(Etol/2)
gam = 50
tic(); Hout3 = HLWB_TwoStep_Accelerated(A,Dgraph,E,F,W,Float64(Etol),Float64(TriTol),lam,"HLWB_new_Tfix.txt",Float64(gam),Float64(1/gam),vareps); HLWB_new_time3 = toc()
hIt3 = Hout3[2]
M = Hout3[1]
D2 = M+M'
HLWB_new_obj3 = sum((A[i,j]-lam)*D2[i,j] for i=1:n-1 for j = i+1:n) + lam*(n*(n-1)/2 - numedges)


# 5. Now the lazy constraints method
tic(); Dtrue, lccbound = FastLPlamCC(A,lam); Exact = toq()

@show lccbound, Exact
@show Dykstra_obj, Dykstra_time, DyIt
@show HLWB_new_obj, HLWB_new_time, hIt
@show HLWB_new_obj2, HLWB_new_time2, hIt2
@show HLWB_new_obj3, HLWB_new_time3, hIt3
