addpath('Community_BGLL_Matlab/')
addpath('GenLouvain-2.1/')
addpath('data')
addpath('Algs')
load Netscience

% Run the louvain algorithm on A
COMTY = cluster_jl_cpp(triu(full(A)),1,1,0,1);
cLou = COMTY.COM{end};

%% Compute the modularity score

Q = modularity(A,cLou)

compute_modularity(cLou,A)

%% A slightly different implementation does some different things with
% randomization, leading to sometimes better objectives

m = nnz(A)/2;
w = sum(A,1);
lam = 1/(2*m);
B = @(i) A(:,i) - lam*w'*w(i);
%B = full(A - w'*w*lam);
limit = 10000;
c = iterated_genlouvain(B,limit);
Q = modularity(A,c)