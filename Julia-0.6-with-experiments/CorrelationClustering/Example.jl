include("Run_CC_exps.jl")


# Choose the name of a graph that is stored in the "../graphs" folder
name = "karate"
Run_DykstraCC(name)

# Run Gurobi Lazy Constraints Method.
# Since we've already generated the instance of correlation clustering
# for karate, we can just load it this time.
generateCC = false
Run_LazyGurobi(name,"experiment_tag",generateCC)
