using MAT
include("software/Gurobi_LR.jl")
name = “karate”

mat = matread("../graphs/"*name*".mat")
A = mat["A"]
n = size(A,1)
# Gurobi subroutines:
# -1=automatic,
# 0=primal simplex,
# 1=dual simplex,
# 2=barrier,
# 3=concurrent,
# 4=deterministic concurrent
# 5=deterministic concurrent simplex.

# Crossover strategy: an integer from -1 to 0. -1 is automatic.
# 0 gets rid of this parameter, which is the fastest.
tic()
X, LR, LPstatus = LeightonRao(A,500,1e-6,2,0)
timeGurobiLR = toc()

##

gam = 20.0
lam = 1/n

include("software/DykstraSC.jl")
maxits = 1000
GapTol = 1e-5
ConTol = 1e-6
tic()
Xdyk, FinalCond, FinalGapDyk, FinalobjDyk, FinalitsDyk, R, LPobjsDyk, duals, primals, gaps,
ConviolationDyk, Bty = Dykstra_LeightonRao_TFA(A,GapTol,ConTol,lam,"DykstraSCoutput",Float64(gam),maxits,10,1e-19)
DykstraTime = toq()

##

gam = 20.0
sig = 1/50
maxits = FinalitsDyk
TriTol = 1e-6
ChangeTol = 1e-5
filename = "BauschkeOutput"
include("software/BauschkeSC.jl")
tic()
Xbau, FinalCon, Finalobj, Finalits, LPobjsBau, QPobjsBau, ConViolationBau, xChange =
BauschkeSC(A,ChangeTol,TriTol,lam,filename,gam,sig,maxits)
BauschkeTime = toq()

