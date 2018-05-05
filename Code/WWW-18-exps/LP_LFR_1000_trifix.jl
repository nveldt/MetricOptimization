
using MAT
include("TriFixSmallMemory.jl")
include("../LamCC_LP.jl")

# A = generate_LFR(n,15,50,.1,20,50)
mat = matread("A1000_exps.mat")
A = mat["A"]

Lams = [logspace(-5,-1,10)' (.15:.1:.95)']
n = size(A,1)

open("TriBound_A_1000.txt","w") do f
        write(f,"Clusterings for A1000_main, Tri Fix Algorithm\n")
end

Dgraph = zeros(Int8,n,n)
for i = 1:n-1
    for j = i+1:n
        if A[i,j] < .1
            Dgraph[j,i] = 1
        end
    end
end

TriTol = .1
Etol = .2

numedges = countnz(A)/2

for i = 1:19
    lam = Lams[i]


    Gams = [1/(lam) 1/(lam) 1/(lam) 1/(lam) 1/(lam) 1/(lam) 1/(lam) 1/lam 1.5/lam 2/lam 3/lam 3/lam 30 20 20 20 20 20 20]

    gam = Gams[i]


    W = lam*ones(n,n)

    for i = 1:n-1
        for j = i+1:n
            if A[i,j] > .1
                W[j,i] = 1-lam
            end
        end
    end
    E = zeros(n,n);
    F = -gam*W;
    P = zeros(n,n)
    Q = zeros(n,n)
    M = TriangleFix(A,Dgraph,E,F,W,P,Q,Float64(gam),Float64(Etol),Float64(TriTol),lam,"TriFix_A1000.txt",gam)
    D = M+M'
    lccbound = sum((A[i,j]-lam)*D[i,j] for i=1:n-1 for j = i+1:n) + lam*(n*(n-1)/2 - numedges)

    println("Done with Lambda = ", lam)
    lamshort = round(lam,5)
    matwrite("D_1000_tri"*"_$lamshort.mat",Dict( "D" => D))
    matwrite("bound_1000_tri"*"_$lamshort.mat",Dict( "lccbound" => lccbound))
    open("TriBound_A_1000.txt","a") do f
        write(f, "$lccbound ")
        write(f, "$lamshort ")
        write(f,"\n")
    end
end
