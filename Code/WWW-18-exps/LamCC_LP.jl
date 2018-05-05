# Functions for solving the LambdaCC LP relaxation on large networks
# up to thousands of nodes

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


# Something weird wrong with this timed version
function TimedFastLPlamCC(A,lam)
    n = size(A,1)
    A = A - diagm(diag(A))

    # Create a model that will use Gurobi to solve
    m = Model(solver=GurobiSolver(OutputFlag = 0))
    @variable(m, x[1:n,1:n])

    # Minimize length of tour
    if lam < .1
    @objective(m, Min, sum((A[i,j]/lam-1)*x[i,j] for i=1:n-1 for j = i+1:n))
    else
    @objective(m, Min, sum((A[i,j]-lam)*x[i,j] for i=1:n-1 for j = i+1:n))
    end

    #tic()
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
                    @constraint(m,x[j,k] - x[i,k] - x[i,j] <= 0)
                end
            end
        end
    end
    #println("Set up first round of constraints: ",toc())

    # Find intial first solution
#    tic()
#    println("Begin an LP solve:")
    solve(m)
#    println("First solve: ", toc())

while true
      D = getvalue(x)'  # x will naturally be upper triangular, but
                        # our check constraints function wants a lower
                        # triangular matrix

     # Store violating tuples in a vector
     violations= Vector{Tuple{Int,Int,Int}}()

    # tic()
     find_violations!(D,violations)

     # Iterate through and add constraints
     numvi = size(violations,1)
    # println("Found ", numvi, " violations: ",toc())

     # Violations (a,b,c) are ordered so that a < b, and (a,b) is the value
     # that needs to be positive in the inequality.
     # The other two may be in a different order, so we need to make sure we are
     # only placing constraints in the upper triangular portion of the matrix.
     # We do so by just calling min and max on pairs of nodes
    # tic()
     for v in violations
         @constraint(m,x[v[1],v[2]] - x[min(v[1],v[3]),max(v[1],v[3])] - x[min(v[2],v[3]),max(v[2],v[3])] <= 0)
     end
     #println("Added in the extra constriants: ", toc())

     if numvi == 0
         break
     end
    # tic()
    # println("Begin an LP solve:")
     solve(m)
    # println("Another solve took: ", toc())
 end

    # Objective score and the Distance matrix
    numedges = countnz(A)/2
    D = getvalue(x)
    lccBound = sum((A[i,j]-lam)*D[i,j] for i=1:n-1 for j = i+1:n) + lam*(n*(n-1)/2 - numedges)
    return D, lccBound
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

    # Minimize length of tour
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

    # Find intial first solution
    solve(m)

while true
      D = getvalue(x)'  # x will naturally be upper triangular, but
                        # our check constraints function wants a lower
                        # triangular matrix

     # Store violating tuples in a vector
     violations= Vector{Tuple{Int,Int,Int}}()
     find_violations!(D,violations)

     # Iterate through and add constraints
     numvi = size(violations,1)

     # Violations (a,b,c) are ordered so that a < b, and (a,b) is the value
     # that needs to be positive in the inequality.
     # The other two may be in a different order, so we need to make sure we are
     # only placing constraints in the upper triangular portion of the matrix.
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

    # Objective score and the Distance matrix
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

    # Minimize length of tour
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

    # Find intial first solution
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
# Given an input matrix A, solve the correlation clustering objective exactly
# on A by setting up an Interger Linear Program and solving it with Gurobi
# This version solves the objective with lazy constraints
function LamCC(A,lam)

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

       tic()
       find_violations!(D,violations)

        # Iterate through and add constraints
       numvi = size(violations,1)
       println("Found ", numvi, " violations: ",toc())


       # Violations (a,b,c) are ordered so that a < b, and (a,b) is the value
       # that needs to be positive in the inequality.
       # The other two may be in a different order, so we need to make sure we are
       # only placing constraints in the upper triangular portion of the matrix.
       # We do so by just calling min and max on pairs of nodes
       #tic()
        for v in violations
            @lazyconstraint(cb,x[v[1],v[2]] - x[min(v[1],v[3]),max(v[1],v[3])] - x[min(v[2],v[3]),max(v[2],v[3])] <= 0)
        end
        #println("Added in the extra constriants: ", toc())

    end

    addlazycallback(m,triangleFix)
    println("Working on solving the ILP currently")
    solve(m)
    println("Done with ILP, working on constraint checking.")


    # Return clustering and objective value
    numedges = countnz(A)/2
    D = getvalue(x)
    Disagreements = sum((A[i,j]-lam)*D[i,j] for i=1:n-1 for j = i+1:n) + lam*(n*(n-1)/2 - numedges)
    return extractClustering(getvalue(x)), Disagreements
end  # end LazyLamCC
