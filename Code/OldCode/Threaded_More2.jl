# At a specific triangle i,j,k, with TriNum indicating which constraint we are
# looking at.
# I follow this edge ordering:
# Edge 1 = ij
# Edge 2 = ik
# Edge 3 = jk
# TriNum tells us which one of these edges is the one that is different from
# the others, e.g. if the constraint is Eij - Eik - Ejk < b, then Eij is the
# "different" variable and TriNum = 1
function OneProjection!(j::Int64,k::Int64,Eij::Float64,Eik::Float64,Ejk::Float64,Dij::Int8,Dik::Int8,Djk::Int8,TriNum::Int64,mu::Float64)
    # Check constraint TriNum for triplet i,j,k where i < j < k
    # There are 3 different triangle inequality constraints for Eij, Eik, Ejk
    # so TriNum could be 1 of 3 numbers
    epsi = 1e-8
    # s = ones(3)
    # s[TriNum] = -1

    if TriNum == 1
        b = -Dij+Dik+Djk
        mu = Eij-Eik-Ejk-b
    elseif TriNum == 2
        b = Dij-Dik+Djk
        mu = -Eij+Eik-Ejk-b
    else
        b = Dij+Dik-Djk
        mu = -Eij-Eik+Ejk-b
    end

    # if TriNum == 1, we're trying to check
    # Eij - Eik - Ejk <= -Dij + Dik + Djk + epsi
    if mu > epsi
        return true
      #push!(local_corr,(TriNum,j,k,mu))
    else
      return false
    end

    # This code is cleaner, but might be a little bit slower. We'll see.
end


function OP2!(j::Int64,k::Int64,Eij::Float64,Eik::Float64,Ejk::Float64,Dij::Int8,Dik::Int8,Djk::Int8,TriNum::Int64,local_corr::Vector{Tuple{Int64,Int64,Int64,Float64}})
    # Check constraint TriNum for triplet i,j,k where i < j < k
    # There are 3 different triangle inequality constraints for Eij, Eik, Ejk
    # so TriNum could be 1 of 3 numbers
    epsi = 1e-8
    # s = ones(3)
    # s[TriNum] = -1

    if TriNum == 1
        b = -Dij+Dik+Djk
        mu = Eij-Eik-Ejk-b
    elseif TriNum == 2
        b = Dij-Dik+Djk
        mu = -Eij+Eik-Ejk-b
    else
        b = Dij+Dik-Djk
        mu = -Eij-Eik+Ejk-b
    end

    # if TriNum == 1, we're trying to check
    # Eij - Eik - Ejk <= -Dij + Dik + Djk + epsi
    if mu > epsi

    end

    # This code is cleaner, but might be a little bit slower. We'll see.
end


function ThreadedTripleLoop!(E::Matrix{Float64},D::Matrix{Int8},corrections::Vector{Vector{Tuple{Int64,Int64,Int64,Float64}}})
    n = size(E,1)
    epsi = 1e-8

    Threads.@threads for i = 1:n-2

       for j = i+1:n-1
       Eij = E[j,i]
       Dij = D[j,i]
        for k = j+1:n
        Dik = D[k,i]
        Djk = D[k,j]
        Eik = E[k,i]
        Ejk = E[k,j]

        # Check triangle i,j,k
        b = Dik + Djk - Dij
        mu = (Eij - Ejk - Eik - b)
        if mu > epsi
          push!(corrections[i],(1,j,k,mu))
          # no further adjustment will be needed for this triplet if this happened
          continue
        end
        # Done checking triangle i,j,k

        # Check triangle i,k,j
        b = -Dik + Djk + Dij
        mu = (-Eij - Ejk + Eik - b)
        if mu > epsi
          push!(corrections[i],(2,j,k,mu))
          continue
        end
        # Done checking triangle i,k,j

        # Check triangle j,k,i
        b = Dik - Djk + Dij
        mu = (-Eij + Ejk - Eik - b)
        if mu > epsi
          push!(corrections[i],(3,j,k,mu))
          continue
        end
        # Done checking triangle j,k,i

    end
    end
    end # end triple for loop

end


function ThreadedTripleLoop2!(E::Matrix{Float64},D::Matrix{Int8},corrections::Vector{Vector{Tuple{Int64,Int64,Int64,Float64}}})
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
    @time ThreadedTripleLoop2!(E,D,corrections)
    TotalViolations = 0
    for i = 1:n
        TotalViolations += length(corrections[i])
    end
    println("Version 2, JobI: $TotalViolations")


    println("Again")

    corrections = Vector{Vector{Tuple{Int64,Int64,Int64,Float64}}}()
    for i = 1:n
        push!(corrections,Vector{Tuple{Int64,Int64,Int64,Float64}}())
    end
    @time ThreadedTripleLoop!(E,D,corrections)
    TotalViolations = 0
    for i = 1:n
        TotalViolations += length(corrections[i])
    end
    println("Version 1: $TotalViolations")

    corrections = Vector{Vector{Tuple{Int64,Int64,Int64,Float64}}}()
    for i = 1:n
        push!(corrections,Vector{Tuple{Int64,Int64,Int64,Float64}}())
    end
    @time ThreadedTripleLoop2!(E,D,corrections)
    TotalViolations = 0
    for i = 1:n
        TotalViolations += length(corrections[i])
    end
    println("Version 2, JobI: $TotalViolations")

end

CheckTriples(parse(Int64,ARGS[1]))
