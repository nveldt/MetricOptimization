TripletMat = Dict{Int64,Float64}()

# We may be able to get away with never needing to go from constraint index t
# to the corresponding triplet (i,j,k), so it's enough to be able to go
# from i,j,k to the index t and use that as a key in a dictionary
function ijk_to_ind(i::Int64,
    j::Int64, k::Int64, n::Int64)
    return (i-1)*n^2 + (j-1)*n + k
end


n = 10

@show haskey(TripletMat,10)
@show get(TripletMat,10,0)

TripletMat[10] = .01

# Put down the default and we're safe
@show get(TripletMat,10,0)
@show TripletMat[10]

its = 1
its += 1
@show its
@show length(TripletMat)

## Scratch work for a single update to a triangle constraint

function TriProject(Xij::Float64,Xik::Float64,Xjk::Float64,i::Int64,j::Int64,k::Int64,n::Int64,delta::Float64)
end
