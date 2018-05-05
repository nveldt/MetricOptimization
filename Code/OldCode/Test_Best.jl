using MAT

include("../Algs/Trifix_v2.jl")
include("Tri_Constraint_Library.jl")
include("Hildreth_Xvar/Experimental_HildrethsTFA.jl")


mat = matread("KarateA.mat")
A = mat["A"]

mat = matread("graphs/SmallNets.mat")
Alesmis = mat["Alesmis"]
Afootball = mat["Afootball"]
Apolbooks = mat["Apolbooks"]
Adolphins = mat["Adolphins"]

A = Alesmis

n = 200
mat = matread("graphs/A$n\_sparse.mat")
A = mat["A"]

n = size(A,1)
D = zeros(Int8,n,n)
for i = 1:n-1
    for j = i+1:n
        if A[i,j] < .1
            D[j,i] = 1
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

gam = 20

# Tolerance for how much we want matrix E to change before we declare convergence
Etol = .01

# Tolerance for how well the triangle constraints hold before declaring convergence
TriTol = .01

##

E = zeros(n,n);
F = -gam*ones(n,n);
P = zeros(n,n)
Q = zeros(n,n)
#
tic()
Finalobj, Finalits, Finaltr, residual = Dykstra_TF(A,D,E,F,W,P,Q,[Etol, Etol/2],TriTol,lam,"short_file","Dykstra_output",Float64(gam),1000)
Basetime = toq()

#

# include("Dykstra_TFA2.jl")
# tic()
# X, Newtr, Newobj, Newits = Dykstra_lamCC_TFA(A,Etol,TriTol,lam,"DykstraLamCCoutput",Float64(gam),1000)
# newtime = toq()
#
# include("Dykstra_TFA3.jl")
# tic()
# X, Newtr3, Newobj3, Newits3 = Dykstra_lamCC_TFA(A,Etol,TriTol,lam,"DykstraLamCCoutput",Float64(gam),1000)
# newtime3 = toq()
#
include("Best_Dykstra_Project/Dykstra_TFA4.jl")
tic()
X, Newtr4, Newobj4, Newits4 = Dykstra_lamCC_TFA(A,Etol,TriTol,lam,"DykstraLamCCoutput",Float64(gam),1000)
newtime4 = toq()

# include("Dykstra_TFA_Dict.jl")
# tic()
# X, NewtrD, NewobjD, NewitsD = Dykstra_lamCC_TFA(A,Etol,TriTol,lam,"DykstraLamCCoutput",Float64(gam),1000)
# newtimeD = toq()
##
# Etol = .05
#
# # Tolerance for how well the triangle constraints hold before declaring convergence
# TriTol = .05
#
# include("Dykstra_TFA_mix.jl")
# tic()
# X, NewtrM, NewobjM, NewitsM = Dykstra_lamCC_TFA(A,Etol,TriTol,lam,"DykstraLamCCoutput",Float64(gam),1000)
# newtimeM = toq()

tic()
Dtrue, lccbound = FastLPlamCC(A,lam)
LPtime = toc()

@show lccbound
# @show Finalobj, Finalits, Finaltr
# @show Newobj3, newtime3, Newits3
@show Newobj4, newtime4, Newits4
# @show NewobjD, newtimeD, NewitsD
# @show NewobjM, newtimeM, NewitsM
