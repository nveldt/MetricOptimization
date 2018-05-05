# Given an input graph, convert it to a correlation clustering problem where
# edge weights are assigned based on the community detection experiments of
# Yubo Wang, Linli Xu, Yucheng Chen, and Hao Wang in the  paper:
# A Scalable Approach for General Correlation Clustering
#
# This is a meaningful way to convert a simple unsigned network into an
# instance of correlation clustering in a more principled way that just
# setting edges to be positive or negative as in cluster editing

include("Tri_Constraint_Library.jl")
include("CodePackage/Tricon_Helper.jl")
include("CodePackage/TFA_Projections_Library.jl")
using MAT

# You can run this on any graph stored in a .mat file that has the
# adjacency matrix stored as a sparse matrix A
function RunTFACC_mediumgraphs(name::String,GapTol::Float64=1e-6,ConTol::Float64=1e-10,delta::Float64=.05,epsi::Float64=.01,gam::Float64=10.0,maxits::Int64=50000)

mat = matread("graphs/"*name*".mat")
A = mat["A"]

A = sparse(A)

n = size(A,1)
C = zeros(n,n)
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

outputstring = "CC_output/Output_"*name*"_TFA.txt"
outputstring2 = "CC_output/Output_"*name*"_TFA_record.txt"

Anew = zeros(n,n)
Dmat = zeros(Float64,n,n)
for i = 1:n-1
    for j = i+1:n
        if C[i,j] > 0
           Anew[j,i] = 1
           Anew[j,i] = 1
       else
           Dmat[j,i] = 1
           Dmat[i,j] = 1
        end
    end
end
Anew = sparse(Anew)

tic()
D, LPobjs, primals, duals, gaps, ConViolation, FinalCon, Finalits,
LPbound, FinalGap, R, Bty = Dykstra_CC_TFA(A,abs.(C),Dmat,GapTol,ConTol,outputstring2,gam,maxits)
DykstraTime = toq()

#LPbound = GeneralCCLP_obj(C,D')

matwrite("CC_output/Output_"*name*"_TFA.mat", Dict( "D" => D, "LPbound" => LPbound,
"DykstraTime" => DykstraTime, "FinalGap" => FinalGap, "FinalCon" => FinalCon, "Finalits" => Finalits,
"gam" => gam, "maxits" => maxits, "GapTol"=> GapTol, "ConTol" => ConTol))

open(outputstring,"w") do f
    write(f,"Output from running the CC Dykstra MCLP method \n")
    write(f, "Graph: "*name*"\n")
    write(f, "Input Parameters: Gamma = $gam, GapTol = $GapTol, ConTol = $ConTol \n")
    write(f, "Output satisfies constraints to within $FinalCon \n")
    write(f, "Relative Duality Gap = $FinalGap \n")
    write(f, "LP Relaxation Obj: $LPbound \n")
    write(f, "Compute Time: $DykstraTime")
    write(f,"\n")
end

end

function RunMany()
Names = ["lesmisA","dolphinsA","lesmisA","adjnounA","BuckGraph","footballA","jazzA","polbooksA","Netscience","celegansmetabolicA","powerA","Acagrqc_conn"]

for i = 1:size(Names,1)
name = Names[i]
RunTFACC_mediumgraphs(name)
end

end
