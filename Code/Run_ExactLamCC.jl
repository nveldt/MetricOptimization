include("Tri_Constraint_Library.jl")
using MAT

# Exactly solve the LambdaCC objective on a very small graph

mat = matread("graphs/A200.mat")
A = mat["A"]

Lams = [logspace(-5,-1,10)' (.15:.1:.95)']
n = size(A,1)
Cs = zeros(n,19)
objs = zeros(19,1)
for i = 1:19
    lam = Lams[i]
    c, Disagreements = LazyLamCC(A,lam)

    Cs[:,i] = c
    objs[i] = Disagreements
end

open("clusterings200.txt","w") do f
    for j = 1:19
        ob = objs[j]
        write(f, "$ob ")
        for t = 1:n
            val = trunc(Int,Cs[t,j])
            write(f,"$val ")
        end
        write(f,"\n")
    end
end
