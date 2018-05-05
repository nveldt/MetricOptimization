# Using the framework of Yubo Wang, Linli Xu, Yucheng Chen, and Hao Wang
#  from paper A Scalable Approach for General Correlation Clustering,
#  construct an instance of correlation clustering by computing similarity
#  scores and then computing from them affinity scores for each pair of nodes

include("Tri_Constraint_Library.jl")
include("CodePackage/TFA_Projections_Library.jl")

using MAT

name = "jazzA"
mat = matread("graphs/"*name*".mat")
A = mat["A"]
n = size(A,1)

C = zeros(n,n)
delta = .05
epsi = .01
for i = 1:n-1
    Ni = findnz(A[i,:])
    for j = i+1:n
        Nj = findnz(A[j,:])
        numerator = length(intersect(Nj[1],Ni[1]))
        denominator = length(union(Ni[1],Nj[1]))
        score = numerator/denominator
        C[i,j] = log((1+score-delta)/(1-score+delta))
        if C[i,j] >= 0.0
            C[i,j] += epsi
        else
            C[i,j] -= epsi
        end
    end
end

C = C+C'
for i = 1:n
    C[i,i] = log((2-delta)/delta)
end
@show minimum(abs.(C))
##

lam = .5
tic()
Dtrue, lccbound = FastLPlamCC(A,lam)
ClusEdit = toc()

# C should have both positive and negative entries
tic()
D,LPbound = LazyGeneralCC(C)
LPtime = toc()

Anew = zeros(n,n)
Dmat = zeros(Float64,n,n)
for i = 1:n-1
    for j = i+1:n
        if C[i,j] > 0
           Anew[j,i] = 1
           Anew[j,i] = 1
       else
           Dmat[j,i] = 1
           Dmat[i,j] = 1
        end
    end
end
Anew = sparse(Anew)

##

tic()
GapTol = 1e-6
ConTol = 1e-10
gam = 25.0
maxits = 100003
Doutput = Dykstra_CC_TFA(A,abs.(C),Dmat,GapTol,ConTol,"DykstraCCoutput",gam,maxits)
DykstraTime = toq()
Ddyk = Doutput[1]

include("CodePackage/Tricon_Helper.jl")
@show ClusEdit, LPtime, DykstraTime

@show GeneralCCLP_obj(C,Ddyk'),GeneralCCLP_obj(C,D)
