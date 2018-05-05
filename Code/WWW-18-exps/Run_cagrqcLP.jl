include("LamCC_LP.jl")
using MAT

path = "/homes/lveldt/GitHubRepos/"
file = "Acagrqc_conn"
mat = matread(path*file*".mat")
A = mat["A"]

# mat = matread("SmallNets.mat")
# A = mat["Alesmis"]

Lams = .95:-.1:.15
lamsize = size(Lams,1)
n = size(A,1)
Cs = zeros(n,lamsize)
objs = zeros(lamsize,1)
for i = 1:lamsize
    lam = Lams[i]
    D, c, lccbound = TimeFastLPlamCC(A,lam)

    Cs[:,i] = c
    objs[i] = lccbound
    println("Done with Lambda = ", lam)

    matwrite("D_cagrqc_conn"*"_$lam.mat",Dict( "D" => D))
    matwrite("c_cagrqc_conn"*"_$lam.mat",Dict( "c" => c, "lccbound" => lccbound))
    # open("D_cagrqc_conn"*"_$lam","w") do f
    #     for jj = 1:n
    #         for ii = 1:n
    #             val = D[ii,jj]
    #         write(f,"$val ")
    #         end
    #         write(f,"\n")
    #     end
    # end

end

open("clusterings_cagrqcconn_pt95.txt","w") do f
    for j = 1:lamsize
        lam = Lams[j]
        ob = objs[j]
        write(f, "$ob ")
        write(f, "$lam ")
        for t = 1:n
            val = trunc(Int,Cs[t,j])
            write(f,"$val ")
        end
        write(f,"\n")
    end
end
