include("Tri_Constraint_Library.jl")
using MAT

# This is code demonstrating how to run experiments such as the ones in 
#  "A correlation clustering based framework for community detection"
#
# Here we load a graph and solve LambdaCC LP relaxation on it to obtain lower
# bounds. Unless the graph is reasonably small, this may take too long to
# converge for very small values of lambda.

# We can run code for ca-GrQc, but this will only work for large values
# For tiny lambda we need to use the Triangle-Fixing Algorithm

# path = "/homes/lveldt/GitHubRepos/"
# file = "Acagrqc_conn"
# mat = matread(path*file*".mat")
# A = mat["A"]

# As an example consider running the procedure on a reasonably small graph
# generated from an LFR benchmark graph generator
mat = matread("graphs/A500.mat")
A = mat["A"]


Lams = [logspace(-5,-1,10)' (.15:.1:.95)']
n = size(A,1)
Cs = zeros(n,19)            # place to store clusterings
objs = zeros(19,1)          # place to store objective scores

# This will take longer and longer as lambda decreases.
# Alternative methods may be needed for tiny values of lambda.

for i = 19:-1:1
    lam = Lams[i]
    D,c, lccbound = TimeFastLPlamCC(A,lam)

    Cs[:,i] = c
    objs[i] = lccbound
    matwrite("c_graph"*"_$lam.mat",Dict( "c" => c, "lccbound" => lccbound))
    println("Done with Lambda = ", lam)

    # Uncomment to store distance matrices if desired
    # matwrite("DistcaGrQc"*"_$lam.mat",Dict( "D" => D))

end
