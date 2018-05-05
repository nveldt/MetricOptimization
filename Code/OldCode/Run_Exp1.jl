using MAT

include("Trifix_v2.jl")
include("Tri_Constraint_Library.jl")

function RunExp1(n::Int64)

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

gam = 50.0
TriTol = .01
DykstraTols = [.1, .05, .025, .01]
BauschkeTols = [.1, .05, .025, .01, .005, .0025]
Dtol = DykstraTols[end]
Btol = BauschkeTols[end]
vareps = Btol/2

Lams = [1/n, .01, .05, .1, .25, .5, .75, .9]
Lams = [.1, .5]
numlams = length(Lams)
maxits = 5000

GraphFile = "fileoutput/N$n\_all\_Lambda"
open(GraphFile, "w") do f
        write(f, "Output both algorithms for graph of size $n\n")
end

for t = 1:numlams

    lam = Lams[t]

    W = lam*ones(n,n)
    for i = 1:n-1
        for j = i+1:n
            if A[i,j] > .1
                W[j,i] = 1-lam
            end
        end
    end

    open(GraphFile, "a") do f
            write(f, "\n\nLambda = $lam\n")
    end

    L = Int64(round(lam*1000))
    Dfile = "fileoutput/N$n/Lam$L\_Dykstra"
    Dfile2 = "fileoutput/N$n/Lam$L\_Dykstra\_short"
    Bfile = "fileoutput/N$n/Lam$L\_Bauschke"
    Bfile2 = "fileoutput/N$n/Lam$L\_Bauschke\_short"

    E = zeros(n,n);
    F = -gam*ones(n,n);
    P = zeros(n,n)
    Q = zeros(n,n)
    tic(); Dout = Dykstra_TF(A,Dgraph,E,F,W,P,Q,DykstraTols,TriTol,lam,Dfile2,Dfile,gam,maxits); Dykstra_time = toc()

    E = zeros(n,n)
    F = -gam*ones(n,n)
    tic(); Bout = Bauschke_TF(A,Dgraph,E,F,W,BauschkeTols,TriTol,lam,Bfile2,Bfile,gam,1/gam,vareps,maxits); Bauschke_time = toc()

    @show Dout
    @show Bout

    Dobj = round(Dout[1],3)
    Dit = Dout[2]
    Dres = round(Dout[4],6)
    Dtritol = round(Dout[3],6)
    Dtime = round(Dykstra_time,2)
    Bobj = round(Bout[1],3)
    Bit = Bout[2]
    Bres = round(Bout[4],6)
    Btritol = round(Bout[3],6)
    Btime = round(Bauschke_time,2)

    open(GraphFile, "a") do f
            # write(f, "Algorithm \t\t Obj \t\t Iters \t\t Residual \t\t Tri-tol \t\t Time \n")
            # write(f, "Dykstra \t\t $Dobj \t\t $Dit \t\t $Dres \t\t $Dtritol \t\t $Dtime\n")
            # write(f, "Bauschke\t\t $Bobj \t\t $Bit \t\t $Bres \t\t $Btritol \t\t $Btime\n")
            write(f, "Algorithm \t Obj \t\t Iters \t Residual \t Tri-tol \t Time \n")
            write(f, "Dykstra \t $Dobj \t $Dit \t $Dres \t $Dtritol \t $Dtime\n")
            write(f, "Bauschke\t $Bobj \t $Bit \t $Bres \t $Btritol \t $Btime\n")
    end

end

end
