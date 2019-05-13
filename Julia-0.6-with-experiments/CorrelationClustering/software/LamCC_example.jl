
include("DykstraCC.jl")

using MAT

name = "KarateA"
mat = matread(homedir()"/data/Graphs/"*name*".mat")
A = mat["A"]
lam = .5
D,W = LamCC_DandW(A,lam)

# Define parameters
GapTol = .01
ConTol = .01
filename = "testoutput"
statusFrequency = 10
gam = 10.0
maxits = 1000

# Run it!
X = Dykstra_CC_TFA(A,W,D,GapTol,ConTol,filename,gam,maxits,statusFrequency)


# Now solve the sparsest cut problem

Xsc = 
