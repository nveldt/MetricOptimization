using MAT
using MatrixNetworks
include("software/Gurobi_SC.jl")
include("software/DykstraSC.jl")
include("software/DykstraSC_Helper.jl")

# Run DykstraSC for a single graph. Store graph-specific stats in a .mat file
# and a text file. Add the basic stats to "outputfile", but don't overwrite it.
function Run_Dykstra_SC(name::String,gam::Float64,outputfile::String,filenote::String,lam::Float64,ConTol::Float64,GapTol::Float64,maxits::Int64)

graphoutput = "output/"*name*"_"*filenote

R = 0
mat = matread("../graphs/"*name*".mat")
A = sparse(mat["A"])
n = size(A,1)
d = sum(A,2)

# A priori approximation guarantee
m = Int64(sum(d)/2)
ap = (1 + 1/(2*gam) + lam*n/(2*gam))
println("Setting up for an approximation of $ap")

tic()
X, FinalCon, FinalGap, Finalobj, Finalits, R, LPobjs, duals, primals, gaps,
Conviolation, Bty = Dykstra_LeightonRao_TFA(A,GapTol,ConTol,lam,graphoutput*"_Dykstra",Float64(gam),maxits,10,1e-19)
DykstraTime = toq()

## Get the a posteriori approximation guarantee for sparsest cut outline in the paper.
tic()
Xbound, bound = LR_QP_postBound(A,X,lam)
ApproxLPTime = toc()

## Lower bound on original LP
LB = Bty - bound/(gam)

ApproxRatio = Finalobj/LB

## This is another way to get an a posteriori guarantee
ApproxGuarantee = round((1 + 1/(2*gam) + lam*n/(2*gam))/(1+R),3)

ratio = min(ApproxGuarantee,ApproxRatio)

maxi = maximum(X)
println("Obj = $Finalobj, Its = $Finalits, Gap = $FinalGap, ConVio = $FinalCon, ApproxRatio = $ratio,\n Time = $DykstraTime")

tim = round(DykstraTime,2)
obj = round(Finalobj,4)
numcon = binomial(n,3)*3
println("$name &  $n & $m & $numcon & $FinalGap & $FinalCon & $obj & $tim & $ratio \\\\")

open(outputfile,"a") do f
    write(f, "\nGraph: "*name*", nodes = $n, edges = $m\ngamma = $gam, lam = $lam \nConTol = $ConTol, GapTol = $GapTol \n")
    write(f, "$name & $n & $m & $numcon & $FinalGap & $FinalCon & $obj & $tim & $ratio \\\\ \n")
end

ComputeTime = DykstraTime;

# Save the output
matwrite(graphoutput*".mat",Dict(
"Finalobj" => Finalobj, "FinalCon" => FinalCon, "Finalits" => Finalits,
"FinalGap" => FinalGap, "X" => X, "ApproxRatio" => ApproxRatio, "LPobjs" => LPobjs, "duals" => duals,
"primals" => primals, "gaps" => gaps, "Conviolation" => Conviolation, "ComputeTime" => ComputeTime, "ApproxLPTime" => ApproxLPTime))

end

function RunSmallerGraphs(filetag::String,gamma::Float64)
Names = ["lesmis","polbooks","karate"]

for i = 1:size(Names,1)
    name = Names[i]
    mat = matread("../graphs/"*name*".mat")
    A = sparse(mat["A"])
    n = size(A,1)
    Run_Dykstra_SC(name,gamma,"Dykstra_SC_all_"*filetag,"Gam5Lam1overN",1/n,1e-8,1e-4,500000)
end

end

function RunMultipleGurobiSC(filetag::String,FullFlag,time_limit::Int64)

Names = ["lesmis","polbooks","karate"]
mkdir(filetag)
outputstring = filetag*"/GurobiSC_small.txt"
for i = 1:size(Names,1)
    name = Names[i]
mat = matread("../graphs/"*name*".mat")
A = mat["A"]

OutputFile = filetag*"/Gurobi_"*name*"_output"

if FullFlag
    tic()
    X, LR = LeightonRao(A,time_limit,1e-8,2,0, OutputFile)
    timeGurobiLR = toc()
else
    tic()
    X, LR = LazyLeightonRao(A,Int64(1e8),1e-8,2,0, OutputFile)
    timeGurobiLR = toc()
end

matwrite(filetag*"/"name*"_Gurobi.mat", Dict( "X" => X, "LR" => LR, "time" => timeGurobiLR))

open(outputstring,"a") do f
    write(f, "Graph: "*name*"\n")
    write(f, "Leighton Rao Relaxation: $LR \n")
    write(f, "Time for Lazy LR method: $timeGurobiLR ")
    write(f,"\n")
end

end

end

# FullFlag = 1 if you want to construct the entire constraint matrix
# FullFlag = 0 if you want to try the Lazy Constraints method, which doesn’t work very well for the Leighton-Rao Sparsest Cut relaxation.

function RunGurobiSC(name::String,filetag::String="",FullFlag::Bool=true,time_limit::Int64=3600,contol::Float64=1e-8,GurobiSolver::Int64 = 2)

outputstring = "output/GurobiLR_"*name*"_"*filetag*".txt"
mat = matread("../graphs/"*name*".mat")
A = mat["A"]

OutputFile = "output/Gurobi_"*name*"_output"

if FullFlag
    tic()
    X, LR = LeightonRao(A,time_limit,contol,GurobiSolver,0, OutputFile)
    timeGurobiLR = toc()
else
    tic()
    X, LR = LazyLeightonRao(A,time_limit,contol,GurobiSolver,0, OutputFile)
    timeGurobiLR = toc()
end

matwrite("output/"name*"_Gurobi.mat", Dict( "X" => X, "LR" => LR, "time" => timeGurobiLR))

open(outputstring,"a") do f
    write(f, "Graph: "*name*"\n")
    write(f, "Leighton Rao Relaxation LP score: $LR \n")
    write(f, "Gurobi Time: $timeGurobiLR ")
    write(f,"\n")
end

end
