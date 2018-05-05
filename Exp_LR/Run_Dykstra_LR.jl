using MAT
using MatrixNetworks
include("Tri_Constraint_Library.jl")
include("TFA_LeightonRao.jl")

# /p/mnt/software/julia-0.6.0/bin/julia

# Run DykstraSC for a single graph. Store graph-specific stats in a .mat file
# and a text file. Add the basic stats to "outputfile", but don't overwrite it.
function Run_Dykstra_LR(name::String,gam::Float64,outputfile::String,filenote::String,lam::Float64,ConTol::Float64,GapTol::Float64,maxits::Int64)

graphoutput = "output/"*name*"_"*filenote

R = 0
mat = matread("../Code/graphs/"*name*".mat")
A = sparse(mat["A"])
n = size(A,1)
d = sum(A,2)

m = Int64(sum(d)/2)
ap = (1 + 1/(2*gam) + lam*n/(2*gam))
println("Setting up for an approximation of $ap")

tic()
X, FinalCon, FinalGap, Finalobj, Finalits, R, LPobjs, duals, primals, gaps,
Conviolation, Bty = Dykstra_LeightonRao_TFA(A,GapTol,ConTol,lam,graphoutput*"_Dykstra",Float64(gam),maxits,10,1e-19)
DykstraTime = toq()

## Try to get an interesting new way to find an approximation guarantee
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

function RunSmallerGraphs(filetag::String)
Names = ["lesmisA","dolphinsA","adjnounA","footballA","jazzA","polbooksA","celegansneuralA","Netscience","celegansmetabolicA"]

for i = 1:size(Names,1)
    name = Names[i]
    mat = matread("../Code/graphs/"*name*".mat")
    A = sparse(mat["A"])
    n = size(A,1)
    Run_Dykstra_LR(name,5.0,"Dykstra_LR_all_"*filetag,"Gam5Lam1overN",1/n,1e-4,1e-4,10000)
end

end



function RunLazyLR_smallgraphs(filetag::String)

Names = ["lesmisA","dolphinsA","adjnounA","footballA","jazzA","polbooksA"]
outputstring = "GurobiLR_all_"*filetag*".txt"
for i = 1:size(Names,1)
    name = Names[i]
mat = matread("../Code/graphs/"*name*".mat")
A = mat["A"]

tic()
X, LR = LazyLeightonRao(A,false,10)
timeGurobiLR = toc()

matwrite("output/"name*"_Gurobi.mat", Dict( "X" => X, "LR" => LR, "time" => timeGurobiLR))

open(outputstring,"a") do f
    write(f, "Graph: "*name*"\n")
    write(f, "Leighton Rao Relaxation: $LR \n")
    write(f, "Time for Lazy LR method: $timeGurobiLR ")
    write(f,"\n")
end

end

end

function RunFullGurobiLR_smallgraphs(filetag::String)

Names = ["lesmisA","dolphinsA","adjnounA","footballA","jazzA","polbooksA"]
outputstring = "GurobiFullLR_all_"*filetag*".txt"
for i = 1:size(Names,1)
    name = Names[i]
mat = matread("../Code/graphs/"*name*".mat")
A = mat["A"]

tic()
X, LR = LeightonRao(A,false,10)
timeGurobiLR = toc()

matwrite("output/"name*"_GurobiFull.mat", Dict( "X" => X, "LR" => LR, "time" => timeGurobiLR))

open(outputstring,"a") do f
    write(f, "Graph: "*name*"\n")
    write(f, "Leighton Rao Relaxation: $LR \n")
    write(f, "Time for Lazy LR method: $timeGurobiLR ")
    write(f,"\n")
end

end

end
