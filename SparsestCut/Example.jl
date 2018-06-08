using MAT
include("software/Gurobi_SC.jl")
include("Run_DykstraSC.jl")
name = "lesmis"

mat = matread("../graphs/"*name*".mat")
A = mat["A"]
n = size(A,1)


# Solve the problem with black-box Gurobi software

# Gurobi subroutines:
# -1=automatic,
# 0=primal simplex,
# 1=dual simplex,
# 2=barrier,
# 3=concurrent,
# 4=deterministic concurrent
# 5=deterministic concurrent simplex.
SolverMethod = 2

# Crossover strategy: an integer from -1 to 0. -1 is automatic.
# 0 means no crossover, which is the fastest.
CrossoverStrategy = 0

timeLimit = 500
ConTol = 1e-6
tic()
X, LPsol, LPstatus = LeightonRao(A,timeLimit,ConTol,SolverMethod,CrossoverStrategy)
timeGurobiLR = toc()

# You can also just call RunGurobiSC, which additionally stores the output
# in the output folder
RunGurobiSC(name)

## The Lazy Gurobi Solver works, but just ends up taking longer.
# Uncomment to note performance of this approach.
# tic()
# Xlazy, LPsolLazy = LazyLeightonRao(A,timeLimit,ConTol,SolverMethod,CrossoverStrategy)
# timeGurobiLR = toc()

## Solve the problem using DykstraSC
include("software/DykstraSC.jl")

gam = 10.0
lam = 1/n
maxits = 1000
GapTol = 1e-5
ConTol = 1e-6
tic()
X_dykstra, ConstraintSatisfaction, FinalGap, DykstraObj, Iterations, R, LPobjs, duals, primals, gaps,
Conviolations, Bty = DykstraSC(A,GapTol,ConTol,lam,"output/DykstraSCoutput",Float64(gam),maxits)
DykstraTime = toq()

println("Gurobi LP solution: $LPsol \nDykstraSC LP score: $DykstraObj ")

## You can plot the convergence of DykstraSC if you wish

using Plots

# Convergence is when the primal and dual objective scores meet,
# and the constraint tolerance is low.
data_labels = ["Primal Obj","Dual Obj", "OPT"]
opt = duals[end]
plot(primals, color=:blue,xlim = (50,size(duals,1)),ylim = (-opt/2,2*opt),xlabel = "Number of Iterations",labels = data_labels, title = "Convergence of DykstraSC")
plot!(duals,color=:green,labels = data_labels)
plot!(opt*ones(size(duals,1),1),color =:red, linestyle =:dash,labels = data_labels)
