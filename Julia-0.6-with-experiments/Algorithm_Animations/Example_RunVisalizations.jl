# Load a graph and create animations of an L2 metric nearness problem
# that is equivalent to a quadratic relaxation of unweighted, complete
# correlation clustering.

using MAT

mat = matread("TinyGraphs.mat")
Kmat = matread("KarateVis.mat")

# Zachary's karate club network
Ak = Kmat["A"]
xyk = Kmat["xy"]

# A2 is a 20 node graph
A2 = mat["A2"]
xy2 = mat["xy2"]

# A1 is a 10 node graph
A1 = mat["A1"]
xy1 = mat["xy1"]

include("MetricNearnessL2_Animations.jl")

## Go all the way through the round before updating the plot
anim1 = RoundMetricLoopAnimate!(10,Ak,xyk)
gif(anim1, "newRound10_Karate_fps2.gif", fps = 2)


## A Hybrid approach: begin by visualizing each step, and then speed up
# and update the plot every n^2 constraints visited

anim2 = HybridMetricLoop!(8,Ak,xyk)
gif(anim2, "Hybrid_Karate_fps5.gif", fps = 5)
gif(anim2, "Hybrid_Karate_fps10.gif", fps = 10)
gif(anim2, "Hybrid_Karate_fps20.gif", fps = 20)

## Go all the way through all steps, updating the plot after each
# visit to a triplet where a nontrivial projection ocurred
anim3 = AllStepsMetricL2!(3,A2,xy2)
gif(anim3, "Allsteps_3rounds_A2_fps10.gif", fps = 10)
