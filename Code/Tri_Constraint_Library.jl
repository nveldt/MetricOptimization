# This library contains functions for solving problems related to correlation
# correlation clustering, including exactly solving LambdaCC, as well as solving
# LP relaxations of LambdaCC.
#
# For exact solutions, we employ the lazy constraints library. We do something
# similar for solving LP relaxations, although there is no actualy
# 'lazyconstraints' feature for this so we just re-solve a new LP with updated
# constraints from scratch each time

using JuMP
using Gurobi

# extractClustering
# Given a 0,1 indicator matrix for node clustering, extract the n x 1
# cluster indicator vector
function extractClustering(x)
    # Assuming the triangle inequality results work, we just have to go through
    # each row of x, and if an element hasn't been placed in a cluster yet, it
    # it starts its own new one and adds everything that is with it.
    n = size(x,1)
    NotClustered = fill(true,n)
    c = zeros(n,1)
    clusnum = 1
    for i = 1:n
        if NotClustered[i]
            for j = i:n
                if x[i,j] < .01 # must be zero, they are clustered together
                    c[j] = clusnum
                    NotClustered[j] = false;
                end
            end
            clusnum +=1
        end
    end
    return c
end

# fiveLP_round
# Rounding procedure for the LP relaxation of the LambdaCC objective.
# Proven to give an answer at most 5 times the optimal solution, as long as
# Lambda >= .5
function fiveLP_round(Dist)
    n = size(Dist,1)

    Vinds = 1:n
    numvinds = n;
    C = zeros(n,1)
    while numvinds > 0
        i = rand(1:numvinds)
        pivot = Vinds[i];
        scores = Dist[pivot,Vinds]

        Tindicator = zeros(numvinds,1)
        for i = 1:numvinds
            if scores[i] <= 2/5
            Tindicator[i] = 1
            end
        end

        T = find(Tindicator)
        Tscores = scores[T]

        average = mean(Tscores)

        if average < 1/5
            NewCinds = [pivot; Vinds[T]]
        else
            NewCinds = pivot
        end

        ci = zeros(n,1)
        ci[NewCinds] = 1
        C = [C ci];
        Vinds = setdiff(Vinds,NewCinds)
        numvinds = size(Vinds)[1]
    end
    C = C[:,2:end]
    m = size(C,2);
    v = (1:m)'
    c = C*v'

    return c
end

# find_violations
# Given a candidate distance matrix D, iterate through all 3-tuples of nodes
# and return the tuples where triangle inequality constraints have been violated.
#
# Output is stored in vector 'violations'
#
# Heuristic speed ups:
#       - Only grab Dij, Dik, Djk once per tuple
#       - if a - b - c > 0, then a>b and a>c. By checking these first
#           we rule out a lot of situations where constraints do not need to be
#           be added, so in all test cases this gave a speed up
#
# Note that we want D to be lower triangular here, though if it is symmetric
# this will work fine as well. We just need to make sure the distance information
# in D is not just stored in the upper triangular portion
function find_violations!(D::Matrix{Float64}, violations::Vector{Tuple{Int,Int,Int}})
  n = size(D,1)

  # We only need this satisfied to within a given tolerance, since the
  # optimization software will only solve it to within a certain tolerance
  # anyways. This can be tweaked if necessary, epsi = 1e-8 yields good results.
  epsi = 1e-8
  @inbounds for i = 1:n-2
       for j = i+1:n-1
          a = D[j,i]
           for k = j+1:n
              b = D[k,i]
              c = D[k,j]
        if a - b > epsi && a - c > epsi && a-b-c > epsi
            push!(violations, (i,j,k))
                # @constraint(m, x[i,j] - x[i,k] - x[j,k] <= 0)
        end
        if b - a > epsi && b - c > epsi && b-a-c > epsi
            push!(violations, (i,k,j))
            # @constraint(m, x[i,k] - x[i,j] - x[j,k] <= 0)
        end

        if c - a > epsi && c-b>epsi && c-a-b > epsi
            push!(violations, (j,k,i))
            # @constraint(m, x[j,k] - x[i,k] - x[i,j] <= 0)
        end
      end
    end
  end
end

function LazyGeneralCC(W)
    n = size(W,1)
    W = W - diagm(diag(W))

    # Create a model that will use Gurobi to solve
    m = Model(solver=GurobiSolver(OutputFlag = 0))
    @variable(m, x[1:n,1:n])

    # Minimize the weight of disagreements (relaxation)
    @objective(m, Min, sum((W[i,j])*x[i,j] for i=1:n-1 for j = i+1:n))

    # Constraints: 0 <= x_ij <= 1 for all node pairs (i,j)
    for i = 1:n-1
        for j = i:n
            @constraint(m,x[i,j] <= 1)
            @constraint(m,x[i,j] >= 0)
        end
    end

    # If we solve with no triangle inequality constraints, we would get the
    # solution D[i,j] = D[j,i] = 1 if A[i,j] = 0 and D[i,j] = D[j,i] = 1
    # otherwise.
    #
    # If we then checked all 3-tuples of nodes, we'd find that the only time
    # triangle constraints were violated is at bad triangles (++- triangles,
    # a tuple that is an unclosed wedge).
    # Hence, from the very beginning we know we must include
    # a triangle inequality constraint at every bad triangle.

    # for i = 1:n
    #     NeighbsI = find(A[:,i])
    #     numNeighbs = size(NeighbsI,1)
    #     for u = 1:numNeighbs-1
    #         j = NeighbsI[u]
    #         for v = u+1:numNeighbs
    #             k = NeighbsI[v]
    #             if A[j,k] == 0
    #                 # Then we have a bad triangle: (i,j), (i,k) \in E
    #                 # but (j,k) is not an edge, so D(j,k) wants to be 1
    #                 #assert(i<j<k)
    #                 @constraint(m,x[min(j,k),max(j,k)] - x[min(i,k),max(i,k)] - x[min(i,j),max(i,j)] <= 0)
    #             end
    #         end
    #     end
    # end

    # Find intial first solution
    solve(m)

    while true
          D = getvalue(x)'  # x will naturally be upper triangular, but
                            # our 'find_violations' function wants a lower
                            # triangular matrix

         # Store violating tuples in a vector
         violations = Vector{Tuple{Int,Int,Int}}()
         find_violations!(D,violations)

         # Iterate through and add constraints
         numvi = size(violations,1)

         # Violations (a,b,c) are ordered so that a < b, and (a,b) is the value
         # that needs to be positive in the inequality:
         # x_ab - x_ac - x_bc <= 0.
         # The other two (b and c) could satisfy either b < c or c <= b.
         # We need to make sure we are only placing constraints in the upper
         # triangular portion of the matrix.
         # We do so by just calling min and max on pairs of nodes
         for v in violations
             #assert(v[1]<v[2])
             @constraint(m,x[v[1],v[2]] - x[min(v[1],v[3]),max(v[1],v[3])] - x[min(v[2],v[3]),max(v[2],v[3])] <= 0)
         end

         if numvi == 0
             break
         end
         solve(m)
     end

    # Return the objective score and the distance matrix
    D = triu(getvalue(x))
    LPbound = 0.0
    for i = 1:n-1
        for j = i+1:n
            if W[i,j] > 0
                LPbound += W[i,j]*D[i,j]
            else
                LPbound += -W[i,j]*(1-D[i,j])
            end
        end
    end
    return D+D', LPbound
end


# FastLPlamCC
# Lazy constraint version for solving the LP relaxation of the LambdaCC
# objective
function FastLPlamCC(A,lam)
    n = size(A,1)
    A = A - diagm(diag(A))

    # Create a model that will use Gurobi to solve
    m = Model(solver=GurobiSolver(OutputFlag = 0))
    @variable(m, x[1:n,1:n])

    # Minimize the weight of disagreements (relaxation)
    @objective(m, Min, sum((A[i,j]-lam)*x[i,j] for i=1:n-1 for j = i+1:n))

    # Constraints: 0 <= x_ij <= 1 for all node pairs (i,j)
    for i = 1:n
        for j = i:n
            @constraint(m,x[i,j] <= 1)
            @constraint(m,x[i,j] >= 0)
        end
    end

    # If we solve with no triangle inequality constraints, we would get the
    # solution D[i,j] = D[j,i] = 1 if A[i,j] = 0 and D[i,j] = D[j,i] = 1
    # otherwise.
    #
    # If we then checked all 3-tuples of nodes, we'd find that the only time
    # triangle constraints were violated is at bad triangles (++- triangles,
    # a tuple that is an unclosed wedge).
    # Hence, from the very beginning we know we must include
    # a triangle inequality constraint at every bad triangle.

    for i = 1:n
        NeighbsI = find(A[:,i])
        numNeighbs = size(NeighbsI,1)
        for u = 1:numNeighbs-1
            j = NeighbsI[u]
            for v = u+1:numNeighbs
                k = NeighbsI[v]
                if A[j,k] == 0
                    # Then we have a bad triangle: (i,j), (i,k) \in E
                    # but (j,k) is not an edge, so D(j,k) wants to be 1
                    #assert(i<j<k)
                    @constraint(m,x[min(j,k),max(j,k)] - x[min(i,k),max(i,k)] - x[min(i,j),max(i,j)] <= 0)
                end
            end
        end
    end

    # Find intial first solution
    solve(m)

    while true
          D = getvalue(x)'  # x will naturally be upper triangular, but
                            # our 'find_violations' function wants a lower
                            # triangular matrix

         # Store violating tuples in a vector
         violations = Vector{Tuple{Int,Int,Int}}()
         find_violations!(D,violations)

         # Iterate through and add constraints
         numvi = size(violations,1)

         # Violations (a,b,c) are ordered so that a < b, and (a,b) is the value
         # that needs to be positive in the inequality:
         # x_ab - x_ac - x_bc <= 0.
         # The other two (b and c) could satisfy either b < c or c <= b.
         # We need to make sure we are only placing constraints in the upper
         # triangular portion of the matrix.
         # We do so by just calling min and max on pairs of nodes
         for v in violations
             #assert(v[1]<v[2])
             @constraint(m,x[v[1],v[2]] - x[min(v[1],v[3]),max(v[1],v[3])] - x[min(v[2],v[3]),max(v[2],v[3])] <= 0)
         end

         if numvi == 0
             break
         end
         solve(m)
     end

    # Return the objective score and the distance matrix
    numedges = countnz(A)/2
    D = getvalue(x)
    lccBound = sum((A[i,j]-lam)*D[i,j] for i=1:n-1 for j = i+1:n) + lam*(n*(n-1)/2 - numedges)
    return D+D', lccBound
end

## Timed version

# TimeFastLPlamCC
# Lazy constraint version for solving the LP relaxation of the LambdaCC
# objective
# This version just prints out statements about how long each step takes
function TimeFastLPlamCC(A,lam)
    n = size(A,1)
    A = A - diagm(diag(A))

    # Create a model that will use Gurobi to solve
    m = Model(solver=GurobiSolver(OutputFlag = 0))
    @variable(m, x[1:n,1:n])

    # Minimize weight of disagreements (relaxation)
    @objective(m, Min, sum((A[i,j]-lam)*x[i,j] for i=1:n-1 for j = i+1:n))

    for i = 1:n
        for j = i:n
            @constraint(m,x[i,j] <= 1)
            @constraint(m,x[i,j] >= 0)
        end
    end

    # If we solve with no triangule inequality constraints, we would get the
    # solution D[i,j] = D[j,i] = 1 if A[i,j] = 0 and D[i,j] = D[j,i] = 1
    # otherwise.
    #
    # If we then checked all 3-tuples of nodes, we'd find that the only time
    # triangle constraints were violated is at bad triangles (a tuple that is
    # an unclosed wedge). Hence, from the very beginning we know we must include
    # a triangle inequality constraint at every bad triangle

    for i = 1:n
        NeighbsI = find(A[:,i])
        numNeighbs = size(NeighbsI,1)
        for u = 1:numNeighbs-1
            j = NeighbsI[u]
            for v = u+1:numNeighbs
                k = NeighbsI[v]
                if A[j,k] == 0
                    # Then we have a bad triangle: (i,j), (i,k) \in E
                    # but (j,k) is not an edge, so D(j,k) wants to be 1
                    #assert(i<j<k)
                    @constraint(m,x[min(j,k),max(j,k)] - x[min(i,k),max(i,k)] - x[min(i,j),max(i,j)] <= 0)
                end
            end
        end
    end

    # Find initial solution
    solve(m)

    while true
          D = getvalue(x)'  # x will naturally be upper triangular, but
                            # our check constraints function wants a lower
                            # triangular matrix

         # Store violating tuples in a vector
         violations= Vector{Tuple{Int,Int,Int}}()

         tic()
         find_violations!(D,violations)

          # Iterate through and add constraints
         numvi = size(violations,1)
         println("Found ", numvi, " violations: ",toq())


         # Violations (a,b,c) are ordered so that a < b, and (a,b) is the value
         # that needs to be positive in the inequality.
         # The other two may be in a different order, so we need to make sure we are
         # only placing constraints in the upper triangular portion of the matrix.
         # We do so by just calling min and max on pairs of nodes
         tic()
          for v in violations
              @constraint(m,x[v[1],v[2]] - x[min(v[1],v[3]),max(v[1],v[3])] - x[min(v[2],v[3]),max(v[2],v[3])] <= 0)
          end
          println("Added in the extra constriants: ", toq())

         if numvi == 0
             break
         end
         tic()
         println("Begin an LP solve:")
         solve(m)
         println("Another solve took: ", toq())
     end

    # Objective score and the Distance matrix
    numedges = countnz(A)/2
    D = getvalue(x)
    lccBound = sum((A[i,j]-lam)*D[i,j] for i=1:n-1 for j = i+1:n) + lam*(n*(n-1)/2 - numedges)
    return D, fiveLP_round(D+D'), lccBound
end


# LazyLamCC
# Given an input matrix A, solve the LambdaCC objective exactly
# on A by setting up an Interger Linear Program and solving it with Gurobi
# This version solves the objective with lazy constraints
function LazyLamCC(A,lam)

    n = size(A,1)
    A = A - diagm(diag(A))

    # Create a model that will use Gurobi to solve
    m = Model(solver=GurobiSolver(OutputFlag = 0))

    @variable(m, x[1:n,1:n], Bin)

    @objective(m, Min, sum((A[i,j]-lam)*x[i,j] for i=1:n-1 for j = i+1:n))

    function triangleFix(cb)
        println("Inside the triangleFix callback")

        D = getvalue(x)'  # x will naturally be upper triangular, but
                          # our check constraints function wants a lower
                          # triangular matrix

       # Store violating tuples in a vector
       violations= Vector{Tuple{Int,Int,Int}}()

       #tic()
       find_violations!(D,violations)

        # Iterate through and add constraints
       numvi = size(violations,1)
       #println("Found ", numvi, " violations: ",toc())

        for v in violations
            @lazyconstraint(cb,x[v[1],v[2]] - x[min(v[1],v[3]),max(v[1],v[3])] - x[min(v[2],v[3]),max(v[2],v[3])] <= 0)
        end

    end

    addlazycallback(m,triangleFix)
    #println("Working on solving the ILP currently")
    solve(m)
    #println("Done with ILP, working on constraint checking.")


    # Return clustering and objective value
    numedges = countnz(A)/2
    D = getvalue(x)
    Disagreements = sum((A[i,j]-lam)*D[i,j] for i=1:n-1 for j = i+1:n) + lam*(n*(n-1)/2 - numedges)
    return extractClustering(getvalue(x)), Disagreements
end  # end LazyLamCC


# WeakLamCCbound
#
# This is another way to get a lower bound on the lambdaCC objective. It comes
# from the LP relaxation used by Ailon, Charikar, and Newman in proving that
# the Pivot algoritm gives a three approximation.
#
# The LP is:
#
#       min \sum_ij x_ij
#       such that  1 <= x_ij + x_ik + x_jk for every bad triangle t = (i,j,k)
#                  0 <= x_ij <= 1 for every edge (i,j)
#
# The idea is that at every bad triangle, we will have to make a mistake at at
# least one of the edges.
#
# If you replace linear constraints with binary constraints on x_ij, you get a
# tighter lower bound on the optimal LamCC objective, but this problem is of
# course much harder to solve.
#
# Commentary: note that for large values of lambda (think .5 or greater) this
# seems to still be a pretty good lower bound on the optimal objective. However,
# it gives a very loose bound for small lambda.
function WeakLamCCbound(A,lam,ExactFlag)

    n = size(A,1)
    A = A - diagm(diag(A))

    # Create a model that will use Gurobi to solve
    m = Model(solver=GurobiSolver(OutputFlag = 0))

    if ExactFlag
        @variable(m, x[1:n,1:n],Bin)
    else
        @variable(m, x[1:n,1:n])
    end

    @objective(m, Min, sum(((1-2*lam)*A[i,j]+lam)*x[i,j] for i=1:n-1 for j = i+1:n))

    for i = 1:n
        @constraint(m, x[i,i] == 0)
        for j = (i+1):n
            @constraint(m, x[i,j] == x[j,i])
            @constraint(m, x[i,j] >= 0)
            @constraint(m, x[i,j] <= 1)
        end
    end

    # For every bad triangle, we must make a mistake at least one edge
    for i = 1:n
        NeighbsI = find(A[:,i])
        numNeighbs = size(NeighbsI,1)
        for u = 1:numNeighbs-1
            j = NeighbsI[u]
            for v = u+1:numNeighbs
                k = NeighbsI[v]
                if A[j,k] == 0   # Then we have a bad triangle!
                                 # i neighbors both j and k, but A[j,k] = 0
                    @constraint(m,x[i,j] + x[i,k] + x[j,k] >= 1)
                end
            end
        end
    end

    solve(m)

    return getvalue(x), getobjectivevalue(m)
end # end WeakLamCCbound


# To add:
#
# 1. Correlation clustering exact solver, not lamCC
# 2. Cluster deletion exact solver
# 3. Cluster deletion LP relaxation
# 4. Cluster deletion other LP for different lower bound


# LazyLeightonRao
# This is the Leighton-Rao relaxation, where we solve on a subset of constraints
# and then update when a triangle inequality constraint isn't satisfied
#
# QuadFlag is a boolean variable indicating whether you want to solve the
# quadratic objective norm(gamma*vecA + d) instead of linear objective vecA'*d
function LazyLeightonRao(A,QuadFlag,gamma)
    n = size(A,1)
    A = A - diagm(diag(A))

    # Create a model that will use Gurobi to solve
    m = Model(solver=GurobiSolver(OutputFlag = 0))
    @variable(m, x[1:n,1:n])

    # Minimize the weight of disagreements (relaxation)
    if QuadFlag
        @objective(m, Min, sum( (A[i,j]*gamma + x[i,j])^2 for i=1:n-1 for j = i+1:n ) )
    else
        @objective(m, Min, sum((A[i,j])*x[i,j] for i=1:n-1 for j = i+1:n))
    end

#   Hard to tell whether or not we really need the bounds between 0 and 1.
#   They seem to help convergence, but I believe adding the xij <= 1 constraint
#   is not technically part of the relaxation.
#   Constraints: 0 <= x_ij <= 1 for all node pairs (i,j)

    for i = 1:n
        for j = i:n
            #@constraint(m,x[i,j] <= 1)
            @constraint(m,x[i,j] >= 0)
        end
    end

    # Constraint: \sum_{ij} x_{ij} = k for some value k, we'll use 1
    # but in Luca Trevisan's notes he uses something else because the definition
    # of sparsest cut he is working with is a multiplicative factor different

    @constraint(m,sum(x[i,j] for i=1:n-1 for j = i+1:n) == n)

    # If we solve with no triangle inequality constraints, we would get the
    # solution D[i,j] = D[j,i] = 1 if A[i,j] = 0 and D[i,j] = D[j,i] = 1
    # otherwise.
    #
    # If we then checked all 3-tuples of nodes, we'd find that the only time
    # triangle constraints were violated is at bad triangles (++- triangles,
    # a tuple that is an unclosed wedge).
    # Hence, from the very beginning we know we must include
    # a triangle inequality constraint at every bad triangle.

    for i = 1:n
        NeighbsI = find(A[:,i])
        numNeighbs = size(NeighbsI,1)
        for u = 1:numNeighbs-1
            j = NeighbsI[u]
            for v = u+1:numNeighbs
                k = NeighbsI[v]
                if A[j,k] == 0
                    # Then we have a bad triangle: (i,j), (i,k) \in E
                    # but (j,k) is not an edge, so D(j,k) wants to be 1
                    #assert(i<j<k)
                    @constraint(m,x[min(j,k),max(j,k)] - x[min(i,k),max(i,k)] - x[min(i,j),max(i,j)] <= 0)
                end
            end
        end
    end

    # Find intial first solution
    solve(m)

    while true
          D = getvalue(x)'  # x will naturally be upper triangular, but
                            # our 'find_violations' function wants a lower
                            # triangular matrix

         # Store violating tuples in a vector
         violations = Vector{Tuple{Int,Int,Int}}()
         find_violations!(D,violations)

         # Iterate through and add constraints
         numvi = size(violations,1)

         # Violations (a,b,c) are ordered so that a < b, and (a,b) is the value
         # that needs to be positive in the inequality:
         # x_ab - x_ac - x_bc <= 0.
         # The other two (b and c) could satisfy either b < c or c <= b.
         # We need to make sure we are only placing constraints in the upper
         # triangular portion of the matrix.
         # We do so by just calling min and max on pairs of nodes
         for v in violations
             #assert(v[1]<v[2])
             @constraint(m,x[v[1],v[2]] - x[min(v[1],v[3]),max(v[1],v[3])] - x[min(v[2],v[3]),max(v[2],v[3])] <= 0)
         end

         if numvi == 0
             break
         end

         tic()
         println("Begin an LP solve:")
         solve(m)
         println("Another solve took: ", toq())
     end

    # Return the objective score and the distance matrix
    D = getvalue(x)
    LP = sum((A[i,j])*D[i,j] for i=1:n-1 for j = i+1:n)
    return D+D', LP
end

# LeightonRao
#
# Full memory leighton rao. Don't run for anything too large
function LeightonRao(A)
    n = size(A,1)
    A = A - diagm(diag(A))

    # Create a model that will use Gurobi to solve
    m = Model(solver=GurobiSolver(OutputFlag = 0))
    @variable(m, x[1:n,1:n])

    # Minimize the weight of disagreements (relaxation)

    @objective(m, Min, sum((A[i,j])*x[i,j] for i=1:n-1 for j = i+1:n))


    #   Constraints: 0 <= x_ij <= 1 for all node pairs (i,j)
    for i = 1:n
        for j = i:n

            # You can get away with not including the xij <= 1 constraint.
            # Seems to help in practice though.
            @constraint(m,x[i,j] <= 1)
            @constraint(m,x[i,j] >= 0)
        end
    end

    # Constraint: \sum_{ij} x_{ij} = k for some value k, we'll use 1
    # but in Luca Trevisan's notes he uses something else because the definition
    # of sparsest cut he is working with is a multiplicative factor different
    @constraint(m,sum(x[i,j] for i=1:n-1 for j = i+1:n) == n)

    # Triangle Inequality constraints
    for i = 1:n-2
        for j = i+1:n-1
            for k = j+1:n

                @constraint(m,x[min(j,k),max(j,k)] - x[min(i,k),max(i,k)] - x[min(i,j),max(i,j)] <= 0)
                @constraint(m,-x[min(j,k),max(j,k)] + x[min(i,k),max(i,k)] - x[min(i,j),max(i,j)] <= 0)
                @constraint(m,-x[min(j,k),max(j,k)] - x[min(i,k),max(i,k)] + x[min(i,j),max(i,j)] <= 0)

            end
        end
    end

    solve(m)

    # Return the objective score and the distance matrix
    D = getvalue(x)
    LP = sum((A[i,j])*D[i,j] for i=1:n-1 for j = i+1:n)
    return D+D', LP
end


# FastLPlamCCnobound
# see what happens if we don't put an upper bound on x_ij
function FastLPlamCCnobound(A,lam)
    n = size(A,1)
    A = A - diagm(diag(A))

    # Create a model that will use Gurobi to solve
    m = Model(solver=GurobiSolver(OutputFlag = 0))
    @variable(m, x[1:n,1:n])

    # Minimize the weight of disagreements (relaxation)
    @objective(m, Min, sum((A[i,j]-lam)*x[i,j] for i=1:n-1 for j = i+1:n))

    # Constraints: 0 <= x_ij <= 1 for all node pairs (i,j)
    for i = 1:n
        for j = i:n
            #@constraint(m,x[i,j] <= 1)
            @constraint(m,x[i,j] >= 0)
        end
    end

    # If we solve with no triangle inequality constraints, we would get the
    # solution D[i,j] = D[j,i] = 1 if A[i,j] = 0 and D[i,j] = D[j,i] = 1
    # otherwise.
    #
    # If we then checked all 3-tuples of nodes, we'd find that the only time
    # triangle constraints were violated is at bad triangles (++- triangles,
    # a tuple that is an unclosed wedge).
    # Hence, from the very beginning we know we must include
    # a triangle inequality constraint at every bad triangle.

    for i = 1:n
        NeighbsI = find(A[:,i])
        numNeighbs = size(NeighbsI,1)
        for u = 1:numNeighbs-1
            j = NeighbsI[u]
            for v = u+1:numNeighbs
                k = NeighbsI[v]
                if A[j,k] == 0
                    # Then we have a bad triangle: (i,j), (i,k) \in E
                    # but (j,k) is not an edge, so D(j,k) wants to be 1
                    #assert(i<j<k)
                    @constraint(m,x[min(j,k),max(j,k)] - x[min(i,k),max(i,k)] - x[min(i,j),max(i,j)] <= 0)
                end
            end
        end
    end

    # Find intial first solution
    solve(m)

    while true
          D = getvalue(x)'  # x will naturally be upper triangular, but
                            # our 'find_violations' function wants a lower
                            # triangular matrix

         # Store violating tuples in a vector
         violations = Vector{Tuple{Int,Int,Int}}()
         find_violations!(D,violations)

         # Iterate through and add constraints
         numvi = size(violations,1)

         # Violations (a,b,c) are ordered so that a < b, and (a,b) is the value
         # that needs to be positive in the inequality:
         # x_ab - x_ac - x_bc <= 0.
         # The other two (b and c) could satisfy either b < c or c <= b.
         # We need to make sure we are only placing constraints in the upper
         # triangular portion of the matrix.
         # We do so by just calling min and max on pairs of nodes
         for v in violations
             #assert(v[1]<v[2])
             @constraint(m,x[v[1],v[2]] - x[min(v[1],v[3]),max(v[1],v[3])] - x[min(v[2],v[3]),max(v[2],v[3])] <= 0)
         end

         if numvi == 0
             break
         end
         solve(m)
     end

    # Return the objective score and the distance matrix
    numedges = countnz(A)/2
    D = getvalue(x)
    lccBound = sum((A[i,j]-lam)*D[i,j] for i=1:n-1 for j = i+1:n) + lam*(n*(n-1)/2 - numedges)
    return D+D', lccBound
end


# Cluster Deletion LP relaxation, full memory
# i.e. we include all constraints up front, rather
# than checking them and adding them if there is a violation
#
# We also use more memory than in necessary--I'm using a full n x n matrix
# of variables for simplicity. This could be improved, but I expect that this
# will not be the bottleneck. In theory we only need m = (# of edges)
# variables, rather than a full n^2. For this application the size of the
# constraint set is what the real bottleneck is though
function CD_relax_fullmemory(A::SparseMatrixCSC{Float64,Int64},ExactFlag)

    n = size(A,1)
    A = A - diagm(diag(A))

    # Create a model that will use Gurobi to solve
    m = Model(solver=GurobiSolver(OutputFlag = 0))

    if ExactFlag
        @variable(m, x[1:n,1:n], Bin)
    else
        @variable(m, x[1:n,1:n])
    end

    # Minimize the weight of disagreements (relaxation)
    @objective(m, Min, sum((A[i,j]*x[i,j] for i=1:n-1 for j = i+1:n)))

    # Constraints: 0 <= x_ij <= 1 for all node pairs (i,j)
    for i = 1:n
        for j = i:n
            @constraint(m,x[i,j] <= 1)
            @constraint(m,x[i,j] >= 0)
            @constraint(m,x[i,j] == x[j,i])
            if A[i,j] == 0
                @constraint(m,x[i,j] == 1)
            end
        end
    end

    # We do not require a full set of metric constraints, given that xij = 1 for
    # nonedges (i,j). We only require a constraint if two or three pairs of
    # nodes in a triplet are actually edges.
    #
    # Thus, we iterate through the graph one node at a time, and then iterate
    # through pairs of neighbors of that node.
    for i = 1:n
        NeighbsI = find(A[:,i])
        numNeighbs = size(NeighbsI,1)
        for u = 1:numNeighbs-1
            j = NeighbsI[u]
            for v = u+1:numNeighbs
                k = NeighbsI[v]
                if A[j,k] == 1
                    @constraint(m,x[min(j,k),max(j,k)] - x[min(i,k),max(i,k)] - x[min(i,j),max(i,j)] <= 0)
                else
                    # Then we have a bad triangle: (i,j), (i,k) \in E
                    # but (j,k) is not an edge, so D(j,k) wants to be 1
                    @constraint(m,1 - x[min(i,k),max(i,k)] - x[min(i,j),max(i,j)] <= 0)
                end
            end
        end
    end

    # Compute solution
    println("Done setting up constraints. Computing LP solution")
    solve(m)

    D = getvalue(x)
    cdBound = sum((A[i,j])*D[i,j] for i=1:n-1 for j = i+1:n)
    return D, cdBound
end

# Given a clustering vector c of the graph with adjacency matrix A, compute the
# associated cluster deletion score
function CD_score(A,c)
    n = size(A,1)
    score = 0
    for i = 1:n
        NeighbsI = find(A[:,i])
        numNeighbs = size(NeighbsI,1)
        for j = 1:numNeighbs
            if c[i] != c[NeighbsI[j]]
                score += 1
            end
        end
    end
    return score/2
end

# Given the LP relaxation for cluster deletion, round it to find a disjoint
# set of cliques in the network
function CDLP_round(A,Dist,cdBound)
    n = size(Dist,1)

    Vinds = 1:n
    numvinds = n;
    C = zeros(n,1)
    while numvinds > 0
        i = rand(1:numvinds)
        pivot = Vinds[i];
        scores = Dist[pivot,Vinds]

        Tindicator = zeros(numvinds,1)
        for i = 1:numvinds
            if scores[i] < .5
            Tindicator[i] = 1
            end
        end

        T = find(Tindicator)
        NewCinds = [pivot; Vinds[T]]
        ci = zeros(n,1)
        ci[NewCinds] = 1
        C = [C ci];
        Vinds = setdiff(Vinds,NewCinds)
        numvinds = size(Vinds)[1]
    end
    C = C[:,2:end]
    m = size(C,2);
    v = (1:m)'
    c = C*v'

    score = CD_score(A,c)

    approx = CD_score(A,c)/cdBound

    println("This clustering is within $approx of the optimal cluster deletion solution")
    return c, score, approx

end

# A way to get an a posteriori bound on the LP score returned by solving the LR-QP
function LR_QP_postBound(A::SparseMatrixCSC{Float64,Int64},Xhat::Matrix{Float64},lam::Float64)
    n = size(A,1)
    A = A - diagm(diag(A))

    # Create a model that will use Gurobi to solve
    m = Model(solver=GurobiSolver(OutputFlag = 0))
    @variable(m, x[1:n,1:n])

    @objective(m, Max, sum((A[i,j]+lam*(1-A[i,j]))*Xhat[j,i]*x[j,i] for i=1:n-1 for j = i+1:n))


    #   Constraints: 0 <= x_ij <= 1 for all node pairs (i,j)
    for i = 1:n
        for j = i:n
            @constraint(m,x[i,j] >= 0)
            @constraint(m,x[i,j] <= n/(n-1))
            @constraint(m,x[i,j]== x[j,i])
        end
    end

    # Constraint: \sum_{ij} x_{ij} = k for some value k, we'll use 1
    # but in Luca Trevisan's notes he uses something else because the definition
    # of sparsest cut he is working with is a multiplicative factor different
    @constraint(m,sum(x[i,j] for i=1:n-1 for j = i+1:n) == n)

    @constraint(m,sum(x[i,j]*A[i,j] for i=1:n-1 for j = i+1:n) <= sum(Xhat[j,i]*A[i,j] for i=1:n-1 for j = i+1:n))

    solve(m)

    # Return the objective score and the distance matrix
    X = getvalue(x)
    bound = sum((A[i,j]+lam*(1-A[i,j]))*Xhat[j,i]*X[j,i] for i=1:n-1 for j = i+1:n)
    return X, bound
end


# A way to get an a posteriori bound on the LP score returned by solving the
# quadratic programming perturbation of the correlation clustering LP relaxation
#
# This isn't as successful as the Leighton-Rao way to get an a posteriori bound
function LamCC_QP_postBound(A::SparseMatrixCSC{Float64,Int64},Xhat::Matrix{Float64},lam::Float64)
    n = size(A,1)
    A = A - diagm(diag(A))
    numedges = countnz(A)/2
    D = zeros(Float64,n,n)
    W = (1-lam)*ones(n,n)
    for i = 1:n-1
        for j = i+1:n
            if A[i,j] < .1
                D[j,i] = 1
                W[j,i] = lam
            end
        end
    end
    Ehat = Xhat - D
    Fhat = abs.(Ehat)

    # Create a model that will use Gurobi to solve
    m = Model(solver=GurobiSolver(OutputFlag = 0))
    @variable(m, e[1:n,1:n])
    @variable(m, f[1:n,1:n])

    # Minimize the weight of disagreements (relaxation)
    @objective(m, Max, sum( Ehat[j,i]*W[j,i]*e[j,i] + Fhat[j,i]*W[j,i]+f[j,i] for i=1:n-1 for j = i+1:n))


    #   Constraints: 0 <= x_ij <= 1 for all node pairs (i,j)
    for i = 1:n-1
        for j = i+1:n
            @constraint(m,e[i,j]<= f[i,j])
            @constraint(m, -e[i,j]<= f[i,j])
            @constraint(m,e[i,j]== e[j,i])
            @constraint(m,f[i,j]== f[j,i])
            @constraint(m,f[j,i] <= 1)
        end
    end

    # Output must have equal to or better LPCC objective score than Xhat
    @constraint(m,sum(W[j,i]*f[j,i] for i=1:n-1 for j = i+1:n) <= sum(W[j,i]*Fhat[j,i] for i=1:n-1 for j = i+1:n) )

    solve(m)

    # Return the objective score and the distance matrix
    F = getvalue(f)
    E = getvalue(e)
    bound = sum( Ehat[j,i]*W[j,i]*E[j,i] + Fhat[j,i]*W[j,i]+F[j,i] for i=1:n-1 for j = i+1:n)
    return F, bound
end
