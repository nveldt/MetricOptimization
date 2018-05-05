# Get as many Leighton Rao relaxations as possible for small graphs
# This is incredibly expensive, but can be done using the lazy leighton rao method.
# I will then compare how quickly we can get results using projection methods

include("Tri_Constraint_Library.jl")

using MAT

# You can run this on any graph stored in a .mat file that has the
# adjacency matrix stored as a sparse matrix A

function RunLazyLR_smallgraphs(name::String)

mat = matread("graphs/"*name*".mat")
A = mat["A"]
@show size(A)

outputstring = "LR_output/Output_"*name*".txt"
@show outputstring

tic()
XLR, LR = LazyLeightonRao(A,false,10)
timeLazyLR = toc()

matwrite("LR_output/Output_"*name*".mat", Dict( "XLR" => XLR, "LR" => LR, "time" => timeLazyLR))

open(outputstring,"w") do f
    write(f, "Graph: "*name*"\n")
    write(f, "Leighton Rao Relaxation: $LR \n")
    write(f, "Time for Lazy LR method: $timeLazyLR ")
    write(f,"\n")
end

end
