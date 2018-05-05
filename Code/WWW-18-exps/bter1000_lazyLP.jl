include("LamCC_LP.jl")
using MAT


# load a matrix
mat = matread("A_bter_1000_exp.mat")
A = mat["A"]

Lams = [logspace(-5,-1,10); collect(.15:.1:.95)]
lamsize = size(Lams,1)
n = size(A,1)


open("btr1000_lazyLP_bounds.txt","w") do f
    write(f,"LP bound for LambdaCC for 1000-node bter network using lazy constraints method.\n")
end


for i = lamsize:-1:1
    lam = Lams[i]
    D, c, lccbound = TimeFastLPlamCC(A,lam)

    println("Done with Lambda = ", lam)


    open("btr1000_lazyLP_bounds.txt","a") do f
        write(f,"$lam \t $lccbound\n")
    end

end
