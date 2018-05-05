# There are a number of things I need to understand in order to get a paralle
# Dykstra's method working. I'll play around with those pieces here.


# corrections is a vector of tuples storing details for a shift in variables
# because of a constraint violation
# corrections = Vector{Tuple{Int64,Int64,Int64,Float64}}()
#
# # I actually need an entire array of such tuples, n of them, where n is a
# # variable
#
# n = 10;
#
# L = Vector{Vector{Tuple{Int64,Int64,Int64,Float64}}}()
#
#
# for i = 1:n
#     push!(L,Vector{Tuple{Int64,Int64,Int64,Float64}}())
# end
#
#
# # now just push a single thing into one of the parts of L
# push!(L[1],(1,2,3,1.0))
#
#
function test_single_thread(m, n)
    x = zeros(m)
    for ele = 1:m
        for i = 1:n
            x[ele] += rand()
        end
    end
    return x
end

function test_multi_thread(m, n)
    x = zeros(m)
    Threads.@threads for ele = 1:m
        for i = 1:n
            x[ele] += rand()
        end
    end
    return x
end

n = 100000

# I want to experiment with different numbers of threads

a = Vector{Vector{Tuple{Int64,Float64}}}()

for i = 1:n
    push!(a,Vector{Tuple{Int64,Float64}}())
end


readme = rand(n)

tic()
Threads.@threads for i = 1:n
    # Thread grabs a part of a to write to
    locala = a[i]
    id = Threads.threadid()

    # all threads read from 'readme' but never change it
    push!(a[i],(id, readme[i]))
end
toc()
