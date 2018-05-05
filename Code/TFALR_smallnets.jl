using MAT
using MatrixNetworks
include("Tri_Constraint_Library.jl")
outputstring = "LR_output/TFA_LR_blah.txt"

Names = ["KarateA","lesmisA","dolphinsA","lesmisA","adjnounA","BuckGraph","footballA","jazzA","polbooksA","Netscience","celegansmetabolicA"]
Names = ["KarateA"]
size(Names)
R = 0
TriTol = 1e-5
GapTol = 1e-6
gam = Float64(5)
maxits = 10000
X = 0

for i = 1:size(Names,1)
name = Names[i]

mat = matread("graphs/"*name*".mat")
A = sparse(mat["A"])

n = size(A,1)
d = sum(A,2)
D = zeros(n,n)
for i = 1:n
    D[i,i] = d[i]
end
L = D-A
E = eigvals(L)
lam2 = E[2]

lami = 1/n
m = Int64(sum(d)/2)
include("CodePackage/TFA_LeightonRao.jl")

##
tic()
X, FinalCon, FinalGap, Finalobj, Finalits, R, LPobjs,
duals, primals, gaps, ConViolation, Bty = Dykstra_LeightonRao_TFA(A,GapTol,TriTol,lami,"LR2",gam,maxits,10,1e-19)
DykstraTime = toq()


## Try to get an interesting new way to find an approximation guarantee
tic()
Xbound2, bound2 = LR_QP_postBound(A,X,lami)
toc()

LB2 = Bty - bound2/(gam)

@show LB2, Finalobj

ApproxRatio = Finalobj/LB2
##
# As long as n \geq 4 and the optimal scaled sparsest cut puts at least 2 nodes
# in each of the two clusters it forms, we have the following lower bound:
ratio = round((1+ 1/(2*gam) + lami*n/(2*gam))/(1+R),10)
@show ratio

println("Obj = $Finalobj, Its = $Finalits, Gap = $FinalGap, ConVio = $FinalCon, ApproxRatio = $ApproxRatio, OtherRatio = $ratio, Time = $DykstraTime")

tim = round(DykstraTime,2)
obj = round(Finalobj,4)

ratio = min(ratio,ApproxRatio)

println("$name &  $n & $m &  &  & & $obj & $tim & $ratio \\\\")


open(outputstring,"a") do f
    #write(f, "Graph: "*name*", nodes = $n, edges = $m, gamma = $gam \n")
    write(f, "$name &  $n & $m &  &  & & $obj & $tim & $ratio \\\\ \n")
end


end
