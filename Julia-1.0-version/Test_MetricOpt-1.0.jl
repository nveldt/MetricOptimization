using MAT


mat = matread("../Julia-0.6-with-experiments/graphs/polbooks.mat")
A = mat["A"]

include("MetricOptimization-1.0.jl")

# Define a few input parameters
GapTol = .01
ConTol = .01
filename = "testoutput"
statusFrequency = 10
gam = 10.0
maxits = 1000
lam = .2

# Solve a relaxed Lambda-Correlation Clustering Problem
Xlamcc = Dykstra_lamCC_TFA(A,lam, GapTol,ConTol,gam,maxits,statusFrequency)

X = Xlamcc[1]              # The matrix of relaxed distances between nodes
FinalCon = Xlamcc[7]       # constraint tolerance satisfied by final output
FinalGap = Xlamcc[10]      # Final duality gap
Finalobj = Xlamcc[9]       # Final linear program objective score


# Solve the sparsest cut relaxation
lam = .1 # Caveat: the parameter lambda means different things for these
         # two different relaxations.
         # See Veldt, Gleich, Wirth WWW 2018 for details on LambdaCC
         # See Veldt, Gleich, Wirth, Saunderson SIMODS 2019 for sparsest cut relaxation details

Xsc = DykstraSC(A,GapTol,ConTol,lam,filename,gam,maxits,statusFrequency)

X = Xsc[1]
FinalCon = Xsc[2]
FinalGap = Xsc[3]
Finalobj = Xsc[4]
