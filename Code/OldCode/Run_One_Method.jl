using MAT
include("Trifix_Library.jl")
include("Tri_Constraint_Library.jl")

function tryDykstra(n::Int64, lam::Float64, Etol::Float64, TriTol::Float64, gam::Float64)

if n == 1
    # Super easy test case
    mat = matread("KarateA.mat")
    A = mat["A"]
end

# Slightly Larger test cases (be sure to make unweighted)
mat = matread("WWW-18-exps/SmallNets.mat")
Alesmis = mat["Alesmis"]
Afootball = mat["Afootball"]
Apolbooks = mat["Apolbooks"]
Adolphins = mat["Adolphins"]

if n == 2
A = spones(Alesmis)
end

if n == 3
    A = Adolphins
end


if n > 10
    mat = matread("graphs/A$n\_sparse.mat")
    A = mat["A"]
end


n = size(A,1)
Dgraph = zeros(Int8,n,n)
for i = 1:n-1
    for j = i+1:n
        if A[i,j] < .1
            Dgraph[j,i] = 1
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

E = zeros(n,n);
F = -gam*ones(n,n);

# correction variables for the E <= F and -E <= F constraints
P = zeros(n,n)
Q = zeros(n,n)

tic()
M = Dykstra_TriangleFix(A,Dgraph,E,F,W,P,Q,Float64(Etol),Float64(TriTol),lam,"Dykstra_Tfix.txt",Float64(gam))
Dykstra_time = toc()

D = M+M'

# Check the bound
numedges = countnz(A)/2
Dykstra_obj = sum((A[i,j]-lam)*D[i,j] for i=1:n-1 for j = i+1:n) + lam*(n*(n-1)/2 - numedges)

@show Dykstra_obj, Dykstra_time
return Dykstra_obj, Dykstra_time

end


function tryHLWB(n::Int64, lam::Float64, Etol::Float64, TriTol::Float64, gam::Float64)

if n == 1
    # Super easy test case
    mat = matread("KarateA.mat")
    A = mat["A"]
end

# Slightly Larger test cases (be sure to make unweighted)
mat = matread("WWW-18-exps/SmallNets.mat")
Alesmis = mat["Alesmis"]
Afootball = mat["Afootball"]
Apolbooks = mat["Apolbooks"]
Adolphins = mat["Adolphins"]

if n == 2
A = spones(Alesmis)
end

if n == 3
    A = Adolphins
end


if n > 10
    mat = matread("graphs/A$n\_sparse.mat")
    A = mat["A"]
end


n = size(A,1)
Dgraph = zeros(Int8,n,n)
for i = 1:n-1
    for j = i+1:n
        if A[i,j] < .1
            Dgraph[j,i] = 1
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

E = zeros(n,n);
F = -gam*ones(n,n);
tic()
M = HLWB_TriangleFix(A,Dgraph,E,F,W,Float64(Etol),Float64(TriTol),lam,"HLWB_Tfix.txt",Float64(gam))
HLWB_time = toc()

D2 = M+M'
# Check the bound
numedges = countnz(A)/2
HLWB_obj = sum((A[i,j]-lam)*D2[i,j] for i=1:n-1 for j = i+1:n) + lam*(n*(n-1)/2 - numedges)

@show HLWB_obj, HLWB_time
return HLWB_obj, HLWB_time



end


function tryLazy(n::Int64, lam::Float64, Etol::Float64, TriTol::Float64, gam::Float64)

if n == 1
    # Super easy test case
    mat = matread("KarateA.mat")
    A = mat["A"]
end

# Slightly Larger test cases (be sure to make unweighted)
mat = matread("WWW-18-exps/SmallNets.mat")
Alesmis = mat["Alesmis"]
Afootball = mat["Afootball"]
Apolbooks = mat["Apolbooks"]
Adolphins = mat["Adolphins"]

if n == 2
A = spones(Alesmis)
end

if n == 3
    A = Adolphins
end


if n > 10
    mat = matread("graphs/A$n\_sparse.mat")
    A = mat["A"]
end


n = size(A,1)
Dgraph = zeros(Int8,n,n)
for i = 1:n-1
    for j = i+1:n
        if A[i,j] < .1
            Dgraph[j,i] = 1
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

# Now the lazy constraints method
tic()
Dtrue, lccbound = FastLPlamCC(A,lam)
lazytime = toq()

@show lccbound, lazytime
return lccbound, lazytime

end
