# Run modularityCC (essentially) by creating a weighted version of correlation clustering

include("DykstraCC.jl")
using MAT

for graph = ["RogetA","Erdos991A","celegansmetabolicA","Harvard500A","SmaGriA","polblogsA","emailA","Vassar85"]
    mat = matread("../graphs/"*graph*".mat")
    A = mat["A"]
    d = sum(A,1)'
    n = size(A,1)

    # The absolute value of the modularity matrix is the weight matrix for CC
    volA = countnz(A)
    W = abs.(A - d*d'/(volA))

    Dmat = zeros(Float64,n,n)
    I = Vector{Int64}()
    J = Vector{Int64}()
    for i = 1:n-1
        for j = i+1:n
            if A[i,j] == 0
               Dmat[j,i] = 1
               Dmat[i,j] = 1
            end
        end
    end

    GapTol = 1e-4
    ConTol = 0.01
    gam = 2.0
    maxits = 1000000
    output = "output/DykstraCC_outputs_for_Gam_"*string(gam)*graph*".txt"

    tic()
    D, LPobjs, primals, duals, gaps, ConViolation, FinalCon, Finalits,
    LPbound, FinalGap, R, Bty, next_corrections = Dykstra_CC_TFA(A,W,Dmat,GapTol,ConTol,output,gam,maxits,10)
    DykstraTime = toq()

    numnext = length(next_corrections)
    outputstring = "output/"*graph*"_"*"Gam_"*string(gam)*"_outputfile.txt"

    Apriori = (1+1/gam)
    Aposteriori = (1+1/gam)/(1+R)

    matwrite("output/"*graph*"_"*"Gam_"*string(gam)*"_DykstraModularity_output.mat", Dict( "D" => D, "LPbound" => LPbound,
    "DykstraTime" => DykstraTime, "FinalGap" => FinalGap, "FinalCon" => FinalCon, "Finalits" => Finalits,
    "gam" => gam, "maxits" => maxits, "GapTol"=> GapTol, "ConTol" => ConTol, "R" => R, "Bty" => Bty,
    "Apriori" => Apriori, "Aposteriori" => Aposteriori))

    open(outputstring,"w") do f
        write(f,"Output from running the CC Dykstra MCLP method \n")
        write(f, "Graph: "*graph*"\n")
        write(f, "Input Parameters: Gamma = $gam, GapTol = $GapTol, ConTol = $ConTol \n")
        write(f, "Output satisfies constraints to within $FinalCon \n")
        write(f, "Relative Duality Gap = $FinalGap \n")
        write(f, "LP Relaxation Obj: $LPbound \n")
        write(f, "Compute Time: $DykstraTime")
        write(f, "A priori approx ratio = $Apriori\n")
        write(f, "A posteriori approx ratio = $Aposteriori\n")
        write(f, "next_corrections stores $numnext nonzero dual variables. \n")
        write(f,"\n")
    end

end
