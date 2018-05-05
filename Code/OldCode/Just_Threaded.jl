function SerialTripleLoop!(E::Matrix{Float64},D::Matrix{Int8},corrections::Vector{Vector{Tuple{Int64,Int64,Int64,Float64}}})
    n = size(E,1)
    epsi = 1e-8

    for i = 1:n-2

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

function TrueSerialTripleLoop!(E::Matrix{Float64},D::Matrix{Int8},corrections::Vector{Tuple{Int64,Int64,Int64,Float64}})
    n = size(E,1)
    epsi = 1e-8

    for i = 1:n-2

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
          push!(corrections,(1,j,k,mu))
          # no further adjustment will be needed for this triplet if this happened
          continue
        end
        # Done checking triangle i,j,k

        # Check triangle i,k,j
        b = -Dik + Djk + Dij
        mu = (-Eij - Ejk + Eik - b)
        if mu > epsi
          push!(corrections,(2,j,k,mu))
          continue
        end
        # Done checking triangle i,k,j

        # Check triangle j,k,i
        b = Dik - Djk + Dij
        mu = (-Eij + Ejk - Eik - b)
        if mu > epsi
          push!(corrections,(3,j,k,mu))
          continue
        end
        # Done checking triangle j,k,i

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

corrections = Vector{Vector{Tuple{Int64,Int64,Int64,Float64}}}()
for i = 1:n
    push!(corrections,Vector{Tuple{Int64,Int64,Int64,Float64}}())
end
@time ThreadedTripleLoop!(E,D,corrections)
TotalViolations = 0
for i = 1:n
    TotalViolations += length(corrections[i])
end

println("Parallel violations: $TotalViolations")

end

CheckTriples(parse(Int64,ARGS[1]))

println("And do it again for good measure")

CheckTriples(parse(Int64,ARGS[1]))
