using MAT
using MatrixNetworks
include("Tri_Constraint_Library.jl")
include("CodePackage/TFA_Projections_Library.jl")

function Run_Dykstra_lamCC(name::String,gam::Float64,outputfile::String,lam::Float64,ConTol::Float64,GapTol::Float64,maxits::Int64)

R = 0
mat = matread("graphs/"*name*".mat")
A = sparse(mat["A"])
n = size(A,1)

m = Int64(countnz(A)/2)
ap = (1 + 1/(gam))
println("Setting up for an approximation of $ap")

statusFreq = 10
stagTol = 1e-9
tic()
X, LPobjs, primals, duals, gaps, Conviolation, FinalCon, Finalits, Finalobj, FinalGap, R, Bty = Dykstra_lamCC_TFA(A,
GapTol,ConTol,Float64(lam),outputfile*"_algOutput",Float64(gam),maxits,statusFreq,stagTol)
DykstraTime = toq()

## Try to get an interesting new way to find an approximation guarantee
tic()
Fbound, bound = LamCC_QP_postBound(A,X,lam)
ApproxLPTime = toc()

LB = Bty - bound/(gam)

@show bound, LB, Finalobj

ApproxRatio = Finalobj/LB
ApproxGuarantee = round((1 + 1/(gam))/(1+R),3)

ratio = min(ApproxGuarantee,ApproxRatio)

println("Obj = $Finalobj, Its = $Finalits, Gap = $FinalGap, ConVio = $FinalCon, ApproxRatio = $ApproxRatio, OtherGuarantee = $ApproxGuarantee Time = $DykstraTime")

tim = round(DykstraTime,2)
obj = round(Finalobj,4)
println("$name &  $n & $m &  &  & & $obj & $tim & $ratio \\\\")

open(outputfile,"a") do f
    write(f, "Obj = $Finalobj, Its = $Finalits, Gap = $FinalGap, ConVio = $FinalCon, ApproxRatio = $ratio, Time = $DykstraTime\n")
    write(f, "Graph: "*name*", nodes = $n, edges = $m, gamma = $gam lambda = $lam \n")
    write(f, "$name &  $n & $m &  &  & & $obj & $tim & $ratio \\\\ \n")
end

ComputeTime = DykstraTime
ApproxRatio = ratio
# Save the output
matwrite("CC_output/Dykstra_LamCC"*"_"*outputfile*".mat",Dict(
"Finalobj" => Finalobj, "FinalCon" => FinalCon, "Finalits" => Finalits,
"FinalGap" => FinalGap, "X" => X, "ApproxRatio" => ApproxRatio, "LPobjs" => LPobjs,
"duals" => duals,"primals" => primals, "gaps" => gaps,
"Conviolation" => Conviolation, "ComputeTime" => ComputeTime, "ApproxLPTime" => ApproxLPTime))

end
