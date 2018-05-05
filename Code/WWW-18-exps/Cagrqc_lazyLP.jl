include("LamCC_LP.jl")
using MAT

path = "/homes/lveldt/GitHubRepos/"
file = "Acagrqc_conn"
mat = matread(path*file*".mat")
A = mat["A"]

# uncomment this to see how it works
# mat = matread("SmallNets.mat")
# A = mat["Adolphins"]

Lams = [logspace(-5,-1,10); collect(.15:.1:.95)]
lamsize = size(Lams,1)
n = size(A,1)


open("cagrqc_conn_lazyLP_bounds.txt","w") do f
    write(f,"LP bound for LambdaCC for (largest connected component of) cagrqc network using lazy constraints method.\n")
end


for i = lamsize:-1:1
    lam = Lams[i]
    D, c, lccbound = TimeFastLPlamCC(A,lam)

    println("Done with Lambda = ", lam)


    open("cagrqc_conn_lazyLP_bounds.txt","a") do f
        write(f,"$lam \t $lccbound\n")
    end

end
