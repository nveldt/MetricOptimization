# Given an input graph, convert it to a correlation clustering problem where
# edge weights are assigned based on the community detection experiments of
# Yubo Wang, Linli Xu, Yucheng Chen, and Hao Wang in the  paper:
# A Scalable Approach for General Correlation Clustering
#
# This is a meaningful way to convert a simple unsigned network into an
# instance of correlation clustering in a more principled way that just
# setting edges to be positive or negative as in cluster editing

include("software/Gurobi_CC.jl")
include("software/DykstraCC.jl")
using MAT


# Given an input graph A, construct from it an instance of correlation clustering
# in the manner outlined originally by Wang et al.
#
# See details in our paper: https://arxiv.org/abs/1806.01678
#
# Original weighting idea comes:
# Yubo Wang, Linli Xu, Yucheng Chen, and Hao Wang.
# A scalable approach for general correlation clustering.
# In International Conference on Advanced Data Mining and Applications, 2013
#
function ConstructWangCC(A::SparseMatrixCSC{Float64,Int64},name::String,delta::Float64=.05,epsi::Float64=.01)
    n = size(A,1)

    defaultval = log((1-delta)/(1+delta))
    C = defaultval*ones(n,n)
    d = sum(A,2)

    # The adjacency matrix squared gives us nodes within a two hop neighborhood
    # This is useful for computing Jaccard coefficients
    A2 = A^2

    for i = 1:n
        println("Constructing CC problem: Done with node $i")
        # Consider all nodes two hops away
        Ni = findnz(A2[i,:])[1]

        for nj = 1:length(Ni)
            j = Ni[nj]
            if i != j
                # jaccard coefficient:
                numerator = A2[i,j]
                denominator = d[i] + d[j] - numerator

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
    end

    Dmat = zeros(Float64,n,n)

    I = Vector{Int64}()
    J = Vector{Int64}()
    for i = 1:n-1
        for j = i+1:n
            if C[i,j] > 0
               # Anew[j,i] = 1
               # Anew[j,i] = 1
               push!(I,i)
               push!(J,j)
           else
               Dmat[j,i] = 1
               Dmat[i,j] = 1
            end
        end
    end
    Amat = sparse(I,J,vec(ones(length(I),1)),n,n)

    matwrite("CCgraphs/WangCC_"*name*".mat", Dict( "Amat" => Amat, "C" => C,
    "Dmat" => Dmat))

    # Return the adjacency matrix for the positive edges, the weights matrix
    # and the anti-adjacency matrix corresponding to "dissimilarity" scores
    # for the metric nearness formulation
    return Amat, C, Dmat
end

# You can run this on any graph stored in a .mat file that has the
# adjacency matrix stored as a sparse matrix A
#
# Parameters:
#   name = name of the graph, stored as string. Should be in "../graphs" folder
#   generateCC: if true, then the weighted graph has already been generated, and is stored in CCgraphs
#               Otherwise, generate it from scratch.
#
#       GapTol = the tolerance for the relative gap between
#                primal and dual functions. Getting converngece to this tolerance
#               means returning a solution that is within factor (1+GapTol)
#               of the optimal solution
#       ConTol = constraint tolerance
#       delta, epsi = parameters for constructing the weighted CC problem
#       gam = the parameter controlling the relationship between the LP and QP
#       maxits = maximum number of iteration to run
function Run_DykstraCC(name::String,filetag::String="",generateCC::Bool=true,gam::Float64=10.0,GapTol::Float64=1e-4,ConTol::Float64=1e-8,delta::Float64=.05,epsi::Float64=.01,maxits::Int64=500000)

outputstring = "output/DykstraCC_final_"*name*filetag
outputstring2 = "output/DykstraCC_output_"*name*filetag

# Generate the WangCC weighted graph or just load it if possible
if generateCC
    tic()
    mat = matread("../graphs/"*name*".mat")
    A1 = mat["A"]
    A,C,Dmat = ConstructWangCC(A1,name,delta,epsi)
    timeGen = toc()
    println("Took $timeGen seconds to create the instance of correlation clustering")
else
    mat = matread("CCgraphs/WangCC_"*name*".mat")
    A = mat["Amat"]
    C = mat["C"]
    Dmat = mat["Dmat"]
end

W = abs.(C)

tic()
D, LPobjs, primals, duals, gaps, ConViolation, FinalCon, Finalits,
LPbound, FinalGap, R, Bty, next_corrections = Dykstra_CC_TFA(A,W,Dmat,GapTol,ConTol,outputstring2,gam,maxits)
DykstraTime = toq()

# ApproxGuarantee = round((1 + 1/(gam))/(1+R),3)

matwrite("output/DykstraCC_"*name*".mat", Dict( "D" => D, "LPbound" => LPbound,
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

# Run several experiments at once. Default setting should work just fine
function RunManyDykstraCC()
Names = ["lesmis","karate","polbooks"]

for i = 1:size(Names,1)
    name = Names[i]
    Run_Wang_DykstraCC(name)
end

end

# Use the Lazy Constraints Gurobi method to solve the problem
# For explanation of parameters, see previous function, and documentation in the Gurobi_CC.jl (may also be helpful to see Gurobi parameter documentation online: http://www.gurobi.com/documentation/8.0/refman/parameters.html)
function Run_LazyGurobi(name::String,filetag::String="",generateCC::Bool=true,time_limit::Int64=604800,FeasTol::Float64=1e-6,SolverMethod::Int64=2, CrossoverStrategy::Int64 =0,delta::Float64=.05,epsi::Float64 =.01)

summaryOutput = "output/Gurobi_final_"*name*filetag
fullOutput = "output/Gurobi_output_"*name*filetag

# Generate the WangCC weighted graph or just load it if possible
if generateCC
    tic()
    mat = matread("../Code/graphs/"*name*".mat")
    A1 = mat["A"]
    A,C,Dmat = ConstructWangCC(A1,name,delta,epsi)
    timeGen = toc()
    println("Took $timeGen seconds to create the instance of correlation clustering")
else
    mat = matread("CCgraphs/WangCC_"*name*".mat")
    A = mat["Amat"]
    C = mat["C"]
    Dmat = mat["Dmat"]
end

# Call black-box software
tic()
D,LPbound = LazyGeneralCC(A,C,time_limit,FeasTol,SolverMethod,CrossoverStrategy,fullOutput)
LPtime = toc()

# Output results
matwrite("output/Gurobi_output_"*name*".mat", Dict( "D" => D, "LPbound" => LPbound, "LPtime" => LPtime))

open(summaryOutput,"w") do f
    write(f, "Graph: "*name*"\n")
    write(f, "LP Relaxation Obj: $LPbound \n")
    write(f, "Compute Time: $LPtime ")
    write(f,"\n")
end

end
