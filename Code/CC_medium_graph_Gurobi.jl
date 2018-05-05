# Get as many correlation clustering relaxations as possible for small graphs

include("Tri_Constraint_Library.jl")
using MAT

# You can run this on any graph stored in a .mat file that has the
# adjacency matrix stored as a sparse matrix A
function RunLazyCC_mediumgraphs(name::String)

mat = matread("graphs/"*name*".mat")
A = mat["A"]

A = sparse(A)

n = size(A,1)
C = zeros(n,n)
delta = .05
epsi = .01
for i = 1:n-1
    Ni = findnz(A[i,:])
    for j = i+1:n
        Nj = findnz(A[j,:])
        numerator = length(intersect(Nj[1],Ni[1]))
        denominator = length(union(Ni[1],Nj[1]))
        score = numerator/denominator
        C[i,j] = log((1+score-delta)/(1-score+delta))
        if C[i,j] > 0.0
            C[i,j] += epsi
        elseif C[i,j] < 0.0
            C[i,j] -= epsi
        elseif C[i,j] == 0.0 && A[i,j] == 1
            C[i,j] = epsi
        else
            C[i,j] = -epsi
        end
    end
end

C = C+C'
for i = 1:n
    C[i,i] = log((2-delta)/delta)
end

outputstring = "CC_output/Output_"*name*".txt"

tic()
D,LPbound = LazyGeneralCC(C)
LPtime = toc()

matwrite("CC_output/Output_"*name*".mat", Dict( "D" => D, "LPbound" => LPbound, "LPtime" => LPtime))

open(outputstring,"w") do f
    write(f, "Graph: "*name*"\n")
    write(f, "LP Relaxation Obj: $LPbound \n")
    write(f, "Compute Time: $LPtime ")
    write(f,"\n")
end

end

function RunABunch()
Names = ["lesmisA","dolphinsA","lesmisA","adjnounA","BuckGraph","footballA","jazzA","polbooksA","Netscience","celegansmetabolicA","powerA","Acagrqc_conn"]

for i = 1:size(Names,1)
name = Names[i]
RunLazyCC_mediumgraphs(name)
end
end
