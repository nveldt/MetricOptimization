using MAT
using MatrixNetworks
include("Tri_Constraint_Library.jl")
include("CodePackage/TFA_LeightonRao.jl")

# /p/mnt/software/julia-0.6.0/bin/julia

function Run_Dykstra_LR(name::String,gam::Float64,outputfile::String,lam::Float64,TriTol::Float64,GapTol::Float64,maxits::Int64)

R = 0
mat = matread("graphs/"*name*".mat")
A = sparse(mat["A"])
n = size(A,1)
d = sum(A,2)

# D = zeros(n,n)
# for i = 1:n
#     D[i,i] = d[i]
# end
# L = D-A
# E = eigvals(L)
# lam2 = E[2]
#
# lami = lam2/n*2*lamfactor

lami = lam

m = Int64(sum(d)/2)
ap = (1 + 1/(2*gam) + lami*n/(2*gam))
println("Setting up for an approximation of $ap")

tic()
X, FinalCon, FinalGap, Finalobj, Finalits, R, LPobjs, duals, primals, gaps,
Conviolation, Bty = Dykstra_LeightonRao_TFA(A,GapTol,TriTol,lami,outputfile*"_algOutput",Float64(gam),maxits,10,1e-19)
DykstraTime = toq()

## Try to get an interesting new way to find an approximation guarantee
tic()
Xbound2, bound2 = LR_QP_postBound(A,X,lami)
ApproxLPTime = toc()

LB2 = Bty - bound2/(gam)

ApproxRatio = Finalobj/LB2
ApproxGuarantee = round((1 + 1/(2*gam) + lami*n/(2*gam))/(1+R),3)

ratio = min(ApproxGuarantee,ApproxRatio)

maxi = maximum(X)
println("Obj = $Finalobj, Its = $Finalits, Gap = $FinalGap, ConVio = $FinalCon, ApproxRatio = $ratio, Time = $DykstraTime, R = $R, maxEl = $maxi")

tim = round(DykstraTime,2)
obj = round(Finalobj,4)
println("$name &  $n & $m &  &  & & $obj & $tim & $ratio \\\\")

open(outputfile,"a") do f
    write(f, "Graph: "*name*", nodes = $n, edges = $m, gamma = $gam \n")
    write(f, "$name &  $n & $m &  &  & & $obj & $tim & $ratio \\\\ \n")
end

ComputeTime = DykstraTime;

# Save the output
matwrite("LR_output/Dykstra_LR"*"_"*outputfile*".mat",Dict(
"Finalobj" => Finalobj, "FinalCon" => FinalCon, "Finalits" => Finalits,
"FinalGap" => FinalGap, "X" => X, "ApproxRatio" => ApproxRatio, "LPobjs" => LPobjs, "duals" => duals,
"primals" => primals, "gaps" => gaps, "Conviolation" => Conviolation, "ComputeTime" => ComputeTime, "ApproxLPTime" => ApproxLPTime))

end
