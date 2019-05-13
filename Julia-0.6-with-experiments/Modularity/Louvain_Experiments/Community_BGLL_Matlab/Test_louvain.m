A = blkdiag(ones(5),ones(5)); A(1,10) = 1; A(10,1) = 1; A(5,6) = 1; A(6,5) = 1;
A = A-diag(diag(A));
spy(A)
%%
self = 1;
debug = 0;
s = 0;
verbose = 1;
M = triu(A);

spy(M)

%%
[COMTY ending] = cluster_jl(M,s,self,debug,verbose);

cLou = COMTY.COM{end};