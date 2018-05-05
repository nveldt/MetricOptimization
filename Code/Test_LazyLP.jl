include("LazyLibrary.jl")

using MAT

mat = matread("A1000.mat")
A = mat["A"]

mat = matread("SmallNets.mat")
Alesmis = mat["Alesmis"]
Afootball = mat["Afootball"]
Apolbooks = mat["Apolbooks"]
Adolphins = mat["Adolphins"]

#A = Afootball

lam = .5
#
#  tic()
# X2, bound2 = LPrelaxlamCC(A,lam)
# time2 = toc()

tic()
X1, bound1 = LazyRelaxedLamCC(A,lam)
time1 = toc()

# tic()
# cExact, Disagreements = LazyLamCC(A,lam)
# time2 = toc()
#
# @show Disagreements
# @show bound1
