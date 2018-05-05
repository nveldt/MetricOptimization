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

function CleanTripleLoop2!(E::Matrix{Float64},D::Matrix{Int8},corrections::Vector{Vector{Tuple{Int64,Int64,Int64,Float64}}})
    n = size(E,1)

for i = 1:n-2
       for j = i+1:n-1
       Eij = E[j,i]
       Dij = D[j,i]
        for k = j+1:n
        Dik = D[k,i]
        Djk = D[k,j]
        Eik = E[k,i]
        Ejk = E[k,j]
        # Check constraint
        # Eij - Eik - Ejk <= -Dij + Dik + Djk (TriNum == 1)
        OP2!(j,k,Eij,Eik,Ejk,Dij,Dik,Djk,1,corrections[i])

        # Check constraint
        # -Eij + Eik - Ejk <= Dij - Dik + Djk (TriNum == 2)
        OP2!(j,k,Eij,Eik,Ejk,Dij,Dik,Djk,2,corrections[i])
        # Check constraint
        # -Eij - Eik + Ejk <= Dij + Dik - Djk (TriNum == 3)
        OP2!(j,k,Eij,Eik,Ejk,Dij,Dik,Djk,3,corrections[i])

    end
    end
    end # end triple for loop
end

function CleanTripleLoop!(E::Matrix{Float64},D::Matrix{Int8},corrections::Vector{Vector{Tuple{Int64,Int64,Int64,Float64}}})
    n = size(E,1)

for i = 1:n-2

       for j = i+1:n-1
       Eij = E[j,i]
       Dij = D[j,i]
        for k = j+1:n
        Dik = D[k,i]
        Djk = D[k,j]
        Eik = E[k,i]
        Ejk = E[k,j]
        mu = 0.0
        # Check constraint
        # Eij - Eik - Ejk <= -Dij + Dik + Djk (TriNum == 1)
        if OneProjection!(j,k,Eij,Eik,Ejk,Dij,Dik,Djk,1,mu)
            push!(corrections[i],(1,j,k,mu))
            continue
        end

        # Check constraint
        # -Eij + Eik - Ejk <= Dij - Dik + Djk (TriNum == 2)
        if OneProjection!(j,k,Eij,Eik,Ejk,Dij,Dik,Djk,2,mu)
            push!(corrections[i],(2,j,k,mu))
            continue
        end
        # Check constraint
        # -Eij - Eik + Ejk <= Dij + Dik - Djk (TriNum == 3)
        if OneProjection!(j,k,Eij,Eik,Ejk,Dij,Dik,Djk,3,mu)
            push!(corrections[i],(3,j,k,mu))
        end

    end
    end
    end # end triple for loop
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
    n = size(E,1)
    epsi = 1e-8
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
       push!(local_corr,(1,j,k,mu))
       # no further adjustment will be needed for this triplet if this happened
       continue
     end
     # Done checking triangle i,j,k

     # Check triangle i,k,j
     b = -Dik + Djk + Dij
     mu = (-Eij - Ejk + Eik - b)
     if mu > epsi
       push!(local_corr,(2,j,k,mu))
       continue
     end
     # Done checking triangle i,k,j

     # Check triangle j,k,i
     b = Dik - Djk + Dij
     mu = (-Eij + Ejk - Eik - b)
     if mu > epsi
         push!(local_corr,(3,j,k,mu))
     end
     # Done checking triangle j,k,i
 end
 end # end inner double for loop

 end

using MAT
function CheckTriples(n)

mat = matread("KarateA.mat")
A = mat["A"]

mat = matread("graphs/SmallNets.mat")
Alesmis = mat["Alesmis"]
Afootball = mat["Afootball"]
Apolbooks = mat["Apolbooks"]
Adolphins = mat["Adolphins"]


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
println("Version 2, defining extra function: $TotalViolations")

corrections = Vector{Vector{Tuple{Int64,Int64,Int64,Float64}}}()
for i = 1:n
    push!(corrections,Vector{Tuple{Int64,Int64,Int64,Float64}}())
end
@time CleanTripleLoop!(E,D,corrections)
TotalViolations = 0
for i = 1:n
    TotalViolations += length(corrections[i])
end
println("Version 3, cleaner code: $TotalViolations")


end

CheckTriples(parse(Int64,ARGS[1]))

println("And do it again for good measure")

CheckTriples(parse(Int64,ARGS[1]))
