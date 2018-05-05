function ThreadedTripleLoopJ!(E::Matrix{Float64},D::Matrix{Int8},corrections::Vector{Vector{Tuple{Int64,Int64,Int64,Float64}}})
    n = size(E,1)
     Threads.@threads for j = 2:n-1
        JobJ!(E,D,j,corrections[j])
    end
end

function JobJ!(E::Matrix{Float64},D::Matrix{Int8},j::Int64,local_corr::Vector{Tuple{Int64,Int64,Int64,Float64}})
    # This is a task performed by a single thread. It is the inner two loops
    # of the triple for loop
    n = size(E,1); epsi = 1e-8

    for i = 1:j-1

        Xij = E[j,i] + D[j,i]
         for k = j+1:n
         Xik = E[k,i] + D[k,i]; Xjk = E[k,j] + D[k,j]

         mu = Xij - Xik - Xjk
         if mu > epsi
           push!(local_corr,(1,j,k,mu))
           continue
         end

         mu = -Xij + Xik - Xjk
         if mu > epsi
           push!(local_corr,(2,j,k,mu))
           continue
         end

         mu = -Xij - Xik + Xjk
         if mu > epsi
             push!(local_corr,(3,j,k,mu))
         end
     end
     end
 end

function ThreadedTripleLoop!(E::Matrix{Float64},D::Matrix{Int8},corrections::Vector{Vector{Tuple{Int64,Int64,Int64,Float64}}})
    n = size(E,1)
     Threads.@threads for i = 1:n-2
        JobI!(E,D,i,corrections[i])
    end
end

function JobI!(E::Matrix{Float64},D::Matrix{Int8},i::Int64,local_corr::Vector{Tuple{Int64,Int64,Int64,Float64}})
    # This is a task performed by a single thread. It is the inner two loops
    # of the triple for loop
    n = size(E,1); epsi = 1e-8
    for j = i+1:n-1
        Xij = E[j,i] + D[j,i]
         for k = j+1:n
         Xik = E[k,i] + D[k,i]; Xjk = E[k,j] + D[k,j]

         mu = Xij - Xik - Xjk
         if mu > epsi
           push!(local_corr,(1,j,k,mu))
           continue
         end

         mu = -Xij + Xik - Xjk
         if mu > epsi
           push!(local_corr,(2,j,k,mu))
           continue
         end

         mu = -Xij - Xik + Xjk
         if mu > epsi
             push!(local_corr,(3,j,k,mu))
         end
     end
     end

 end

using MAT

function CheckTriples(n)

    mat = matread("KarateA.mat")
    A = mat["A"]

    mat = matread("graphs/SmallNets.mat")
    A = mat["Afootball"]

    mat = matread("graphs/A$n\_sparse.mat")
    A = mat["A"]

    n = size(A,1)
    Dgraph = zeros(Int8,n,n)
    for i = 1:n-1
        for j = i+1:n
            if A[i,j] < .1
                Dgraph[j,i] = 1
            end
        end
    end
    D = Dgraph

    E = zeros(n,n);
    @show Threads.nthreads()

    corrections = Vector{Vector{Tuple{Int64,Int64,Int64,Float64}}}()
    for i = 1:n
        push!(corrections,Vector{Tuple{Int64,Int64,Int64,Float64}}())
    end
    @time ThreadedTripleLoopJ!(E,D,corrections)
    TotalViolations = 0
    for i = 1:n
        TotalViolations += length(corrections[i])
        #@show length(corrections[i])
    end
    println(" JobJ: $TotalViolations")

    corrections = Vector{Vector{Tuple{Int64,Int64,Int64,Float64}}}()
    for i = 1:n
        push!(corrections,Vector{Tuple{Int64,Int64,Int64,Float64}}())
    end
    @time ThreadedTripleLoop!(E,D,corrections)
    TotalViolations = 0
    for i = 1:n
        TotalViolations += length(corrections[i])
        #@show length(corrections[i])
    end
    println("Version JobI: $TotalViolations")


    corrections = Vector{Vector{Tuple{Int64,Int64,Int64,Float64}}}()
    for i = 1:n
        push!(corrections,Vector{Tuple{Int64,Int64,Int64,Float64}}())
    end
    @time ThreadedTripleLoopJ!(E,D,corrections)
    TotalViolations = 0
    for i = 1:n
        TotalViolations += length(corrections[i])
        #@show length(corrections[i])
    end
    println(" JobJ: $TotalViolations")

end

CheckTriples(parse(Int64,ARGS[1]))
