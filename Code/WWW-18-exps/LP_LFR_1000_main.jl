
using MAT
include("TriFixSmallMemory.jl")
include("../LamCC_LP.jl")

# A = generate_LFR(n,15,50,.1,20,50)
mat = matread("A1000_exps.mat")
A = mat["A"]

Lams = [logspace(-5,-1,10)' (.15:.1:.95)']
n = size(A,1)

open("clusteringsA_main_lccbounds.txt","w") do f
        write(f,"Clusterings for A1000_main\n")
end

for i = 1:19
    lam = Lams[i]
    D,c, lccbound = TimeFastLPlamCC(A,lam)

    println("Done with Lambda = ", lam)
    lamshort = round(lam,5)
    matwrite("D_1000_main"*"_$lamshort.mat",Dict( "D" => D))
    matwrite("c_1000_main"*"_$lamshort.mat",Dict( "c" => c, "lccbound" => lccbound))
    open("clusA_main_lccbounds.txt","a") do f
        write(f, "$lccbound ")
        write(f, "$lamshort ")
        for t = 1:n
            val = trunc(Int,c[t])
            write(f,"$val ")
        end
        write(f,"\n")
    end
end
