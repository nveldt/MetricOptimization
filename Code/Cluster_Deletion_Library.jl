using JuMP
using Gurobi

# Cluster Deletion LP relaxation, full memory
# i.e. we include all constraints up front, rather
# than checking them and adding them if there is a violation
#
# This is inefficient and probably shouln't be used.
# There is no need to use a full n x n matrix x instead of a linearized
# vector with an entry for each edge in the graph.
#
# I coded it this way at first for a simple first version, and didn't think
# the bottleneck would be the O(n^2) memory requirement of the variables.
# It may not end up being the bottleneck, but it still needs to a significant
# time cost that is simple to avoid.
function CD_relax_fullmemory(A::SparseMatrixCSC{Float64,Int64},ExactFlag)

    n = size(A,1)
    A = A - diagm(diag(A))

    tic()
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
    time = toq()
    println("Done setting up constraints in $time seconds ")
    solve(m)

    D = getvalue(x)
    cdBound = sum((A[i,j])*D[i,j] for i=1:n-1 for j = i+1:n)
    return D, cdBound
end

# Cluster deletion LP relaxation, version 2
# This version linearizes the x variables so that we only store O(m) variables
# instead of O(n^2).
#
# We still include the entire constraint set up front
function CD_relax_Gurobi(A::SparseMatrixCSC{Float64,Int64})

    n = size(A,1)
    A = A - diagm(diag(A))

    # Always look at lower triangle, where j > i
    I,J,V = findnz(tril(A))

    numedges = length(I)

    # This can now be used to translate from ij to a unique index
    # We could also use a dictionary, which we'll try in the future
    ID = sparse(I,J,collect(1:numedges))

    # Create a model that will use Gurobi to solve
    m = Model(solver=GurobiSolver(OutputFlag = 0))

    @variable(m, x[1:numedges])

    # Minimize the weight of disagreements (relaxation)
    @objective(m, Min, sum((A[I[t],J[t]]*x[t] for t = 1:numedges)))

    # Constraints: 0 <= x_ij <= 1 for all node pairs (i,j)
    for t = 1:numedges
            @constraint(m,x[t] <= 1)
            @constraint(m,x[t] >= 0)
    end

    # Store a list of neighbors for each nodes, since we will frequencly access
    # this and we do not wish to repeatedly call the findnz function.
    # We could avoid storing this, but it shoulnd't be a huge extra memory constraint
    AdjList = Array{Array{Int64}}(n)
    deg = Array{Int64}(n)

    NumConstraints = 0
    for i = 1:n
         AdjList[i] = find(A[:,i])
         deg[i] = length(AdjList[i])
         if deg[i] > 1
             NumConstraints += binomial(deg[i],2)
         end
    end
    println("There are $NumConstraints Constraints for this LP")

    # We do not require a full set of metric constraints, given that xij = 1 for
    # nonedges (i,j). We only require a constraint if two or three pairs of
    # nodes in a triplet are actually edges.
    #
    # Thus, we iterate through the graph one node at a time, and then iterate
    # through pairs of neighbors of that node.
    for i = 1:n
        NeighbsI = AdjList[i]
        numNeighbs = deg[i]

        # Check through pairs of iterates
        for u = 1:numNeighbs-1
            j = NeighbsI[u]
            ij = ID[max(j,i),min(j,i)]
            for v = u+1:numNeighbs      # note that u < v and j < k
                k = NeighbsI[v]
                ik = ID[max(k,i),min(k,i)]
                if A[k,j] == 0
                    @constraint(m,1 - x[ik] - x[ij] <= 0)
                else
                    jk = ID[k,j]
                    @constraint(m,x[jk] - x[ik] - x[ij] <= 0)
                end
            end
        end
    end

    # Compute solution
    println("Done setting up constraints")
    solve(m)

    D = getvalue(x)
    OneMinusD = sparse(I,J,1-D,n,n)
    OneMinusD = OneMinusD + OneMinusD'
    cdBound = sum((A[I[t],J[t]]*D[t] for t = 1:numedges))
    return D, OneMinusD, cdBound

end

# Find violations to the constraint set
function find_CD_violations!(A::SparseMatrixCSC{Float64,Int64},X::Array{Float64},violations::Vector{Tuple{Int64,Int64,Int64}})
    println("Finding violations")
end

# Cluster Deletion LP relaxation, lazy updates version
function Lazy_CD_relax(A::SparseMatrixCSC{Float64,Int64})

    n = size(A,1)
    A = A - diagm(diag(A))

    # Always look at lower triangle, where j > i
    I,J,V = findnz(tril(A))

    numedges = length(I)

    # This can now be used to translate from ij to a unique index
    # We could also use a dictionary, which we'll try in the future
    ID = sparse(I,J,collect(1:numedges))

    # Create a model that will use Gurobi to solve
    m = Model(solver=GurobiSolver(OutputFlag = 0))

    @variable(m, x[1:numedges])

    # Minimize the weight of disagreements (relaxation)
    @objective(m, Min, sum((A[I[t],J[t]]*x[t] for t = 1:numedges)))

    # Constraints: 0 <= x_ij <= 1 for all node pairs (i,j)
    for t = 1:numedges
            @constraint(m,x[t] <= 1)
            @constraint(m,x[t] >= 0)
    end

    # Store a list of neighbors for each nodes, since we will frequencly access
    # this and we do not wish to repeatedly call the findnz function.
    # We could avoid storing this, but it shoulnd't be a huge extra memory constraint
    AdjList = Array{Array{Int64}}(n)
    NumNeigh = Array{Int64}(n)

    for i = 1:n
         AdjList[i] = find(A[:,i])
         NumNeigh[i] = length(AdjList[i])
    end

    # We do not require a full set of metric constraints, given that xij = 1 for
    # nonedges (i,j). We only require a constraint if two or three pairs of
    # nodes in a triplet are actually edges.
    #
    # Thus, we iterate through the graph one node at a time, and then iterate
    # through pairs of neighbors of that node.
    for i = 1:n
        NeighbsI = AdjList[i]
        numNeighbs = NumNeigh[i]

        for u = 1:numNeighbs-1
            j = NeighbsI[u]
            tij = ID[max(j,i),min(j,i)]
            for v = u+1:numNeighbs      # note that u < v and j < k
                k = NeighbsI[v]
                tik = ID[max(k,i),min(k,i)]
                if A[k,j] == 0
                    @constraint(m,1 - x[tik] - x[tij] <= 0)
                end
            end
        end
    end

    # Compute solution
    println("Beginning first solve")
    solve(m)

    while true
          X = getvalue(x)

         # Store violating tuples in a vector
         violations = Vector{Tuple{Int64,Int64,Int64}}()
         find_CD_violations!(A,X,violations)

         # Iterate through and add constraints
         numvi = length(violations)

         # Each violation is stored as (i,j,k) where a is the pivot node
         # and j < k. Index i is not necessarily smaller than j or k.
         for v in violations
             ij = ID[max(j,i),min(j,i)]
             ik = ID[max(k,i),min(k,i)]
             if A[k,j] > 0
                 jk = ID[k,j]
                 @constraint(m,x[jk] - x[ik] - x[ij] <= 0)
             else
                 @constraint(m,1 - x[ik] - x[ij] <= 0)
             end
         end
         if numvi == 0
             break
         end
         solve(m)
     end


    D = getvalue(x)
    cdBound = sum((A[j,i])*D[ID[j,i]] for i=1:n-1 for j = i+1:n)
    return D, I, J, cdBound
end

# Given a clustering vector c of the graph with adjacency matrix A, compute the
# associated cluster deletion score
function CD_score(A::SparseMatrixCSC{Float64,Int64},c::Array{Int64})
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
#
# OneMinusD is a sparse matrix storing the scores 1-xij where xij is the
# "distance" between nodes i and j where (i,j) \in E.
# I store 1-xij instead of xij, because then 1-xuv = 0 then for any non-edge
# (u,v) \notin E.
function CDLP_round_slow(A::SparseMatrixCSC{Float64,Int64},OneMinusD::SparseMatrixCSC{Float64,Int64},cdBound::Float64)
    n = size(A,1)

    Vinds = 1:n
    numvinds = n;
    C = zeros(n,1)

    while numvinds > 0

        # Randomly select an unclustered node
        tic()
        i = rand(1:numvinds)

        # Its global index is "pivot"
        pivot = Vinds[i];

        # Get scores for its neighbors
        scores = OneMinusD[pivot,Vinds]
        toc()

        tic()
        Tindicator = zeros(numvinds,1)
        for i = 1:numvinds
            if scores[i] > .5
            Tindicator[i] = 1
            end
        end
        toc()

        tic()
        T = find(Tindicator)
        NewCinds = [pivot; Vinds[T]]
        ci = zeros(n,1)
        ci[NewCinds] = 1
        C = [C ci];
        Vinds = setdiff(Vinds,NewCinds)
        numvinds = size(Vinds)[1]
        toc()
    end

    tic()
    C = C[:,2:end]
    m = size(C,2);
    v = (1:m)'
    c = C*v'

lastpart= toc()
    score = CD_score(A,c)


    approx = score/cdBound

    println("This clustering is within $approx of the optimal cluster deletion solution")
    return c, score, approx

end


# Given the LP relaxation for cluster deletion, round it to find a disjoint
# set of cliques in the network
#
# OneMinusD is a sparse matrix storing the scores 1-xij where xij is the
# "distance" between nodes i and j where (i,j) \in E.
# I store 1-xij instead of xij, because then 1-xuv = 0 then for any non-edge
# (u,v) \notin E.
function CDLP_round(A::SparseMatrixCSC{Float64,Int64},OneMinusD::SparseMatrixCSC{Float64,Int64},cdBound::Float64)
    n = size(A,1)

    Vinds = 1:n
    numvinds = n;
    C = zeros(n,1)
    while numvinds > 0

        # Randomly select an unclustered node
        tic()
        # Its global index is "pivot"
        pivot = rand(Vinds);

        # Get scores for its neighbors
        scores = OneMinusD[pivot,Vinds]
        toc()

        tic()
        # Tindicator = zeros(numvinds,1)
        # for i = 1:numvinds
        #     if scores[i] > .5
        #     Tindicator[i] = 1
        #     end
        # end
        T= find(scores .> .5)
        toc()

        tic()
        #T = find(Tindicator)
        NewCinds = [pivot; Vinds[T]]
        ci = zeros(n,1)
        ci[NewCinds] = 1
        C = [C ci];
        Vinds = setdiff(Vinds,NewCinds)
        numvinds = size(Vinds)[1]
        toc()
    end

    tic()
    C = C[:,2:end]
    m = size(C,2);
    v = (1:m)'
    c = C*v'

lastpart= toc()
    score = CD_score(A,c)


    approx = score/cdBound

    println("This clustering is within $approx of the optimal cluster deletion solution")
    return c, score, approx

end


# Given the LP relaxation for cluster deletion, round it to find a disjoint
# set of cliques in the network
#
# This function works for dense n x n matrices Dist where the xij = 1
# is explicitly stored
function CDLP_dense_round(A,Dist,cdBound)
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
