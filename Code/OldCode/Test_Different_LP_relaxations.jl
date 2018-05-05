## Test the different LP relaxations
include("Tri_Constraint_Library.jl")

using MAT
mat = matread("graphs/A300.mat")
A = mat["A"]


mat = matread("graphs/SmallNets.mat")
Alesmis = mat["Alesmis"]
Afootball = mat["Afootball"]
Apolbooks = mat["Apolbooks"]
Adolphins = mat["Adolphins"]

A = Apolbooks
n = size(A,1)

lam = .5;   # check dolphins with lam = .1 to see the LP relaxations aren't the same, even if you solve the problem exactly
tic()
D, lccbound1 = FastLPlamCC(A,lam)
toc()
ExactFlag = true;

tic()
D2, bound2 = WeakLamCCbound(A,lam,ExactFlag)
toc()


@show bound2
@show lccbound1
