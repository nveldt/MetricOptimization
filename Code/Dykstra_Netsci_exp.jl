using MAT
using MatrixNetworks

outputstring = "LR_output/TFA_LR_Netscience.txt"

Names = ["Netscience"]

size(Names)
R = 0
TriTol = 1e-5
GapTol = 1e-6
gam = Float64(5)
maxits = 10000

for gam = [1 2 5 10 20]

for i = 1 #1:size(Names,1)

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
lami = lam2/n*20
m = Int64(sum(d)/2)
include("CodePackage/TFA_LeightonRao.jl")

tic()
X, Finaltr, FinalGap, Finalobj, Finalits, R = Dykstra_LeightonRao_TFA(A,GapTol,TriTol,lami,"LR2",Float64(gam),maxits,10,1e-19)
DykstraTime = toq()

ratio = round((1+1/gam)/(1+R),3)
println("Obj = $Finalobj, Its = $Finalits, Gap = $FinalGap, ConVio = $Finaltr, ApproxRatio = $ratio, Time = $DykstraTime")

tim = round(DykstraTime,2)
obj = round(Finalobj,4)
println("$name &  $n & $m &  &  & & $obj & $tim & $ratio \\\\")

open(outputstring,"a") do f
    write(f, "Graph: "*name*", nodes = $n, edges = $m, gamma = $gam \n")
    write(f, "$name &  $n & $m &  &  & & $obj & $tim & $ratio \\\\ \n")
end

end

end
