using MAT
using Plots

include("Cluster_Deletion_Library.jl")

# Run the Gurobi version of the Cluster Deletion LP relaxation
# Save and output stats
function Run_CDrelax_Gurobi(name::String)

    mat = matread("graphs/"*name*".mat")
    A = sparse(mat["A"])

    tic()
    D, OneMinusD, cdBound = CD_relax_Gurobi(A)
    timeCDrelax = toq()
    println("Took $timeCDrelax to compute LP relaxation.")

    tic()
    c, score, approx = CDLP_round(A,OneMinusD,cdBound)
    roundTime = toq()
    println("Took $roundTime seconds to round")
    matwrite("CD_output/CDrelaxGurobi"*"_"*name*".mat",Dict(
    "Dist" => D, "OneMinusDist" => OneMinusD, "compute_time" => timeCDrelax,
    "c" => c,"cd_score" => score, "cd_bound" => cdBound,"round_time" => roundTime))


end
