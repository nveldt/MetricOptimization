using MAT

include("Tri_Constraint_Library.jl")
include("CodePackage/TFA_Projections_Library.jl")
include("CodePackage/TFA_Library_Helper.jl")


mat = matread("graphs/StocksCC_exp.mat")
A = mat["A"]
W = mat["W"]
D,Wnot = LamCC_DandW(A,.1)

gam = 100
GapTol = 1e-1
ConTol = 1e-2
MaxIts = Int64(1500)
statusFreq = 10
stagTol = 1e-9
tic()
X, LPobjs, primals, duals, gaps, ConViolation, Finalits, Finalobj = Dykstra_CC_TFA(A,W,D,
GapTol,ConTol,"DykstraStocksCCoutput",Float64(gam),MaxIts,statusFreq,stagTol)
newtime = toq()
