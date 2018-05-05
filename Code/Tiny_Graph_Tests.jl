include("Tri_Constraint_Library.jl")

n = 3

A = zeros(n,n)

for j = 1:n-1
    A[j,j+1] = 1
    A[j+1,j] = 1
end
# A[n,1]=1
# A[1,n] = 1

A = zeros(n,n)

for i = 1:n
    for j = i+1:n
        if rand(1)[1] > 0
            A[i,j] = 1
            A[j,i] = 1
        end
    end
end

A = sparse(A)
tic()
X, LP = LazyLeightonRao(A,false,75)
time1 = toc()

@show LP
@show maximum(X), countnz(X)/2, n/(n-1), round(maximum(X)*countnz(X),1)/2
