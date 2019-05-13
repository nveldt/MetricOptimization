function [cBest,BestMod,Mods] = many_louvain(A,k,cinit)
% Run the gen-louvain algorithm over and over with different settings
% and take the best outcome
addpath('GenLouvain-2.1/')

n = size(A,1);
if nargin < 3
    cinit = (1:n)';
end

m = nnz(A)/2;
w = sum(A,1);
lam = 1/(2*m);
B = @(i) A(:,i) - lam*w'*w(i);
%B = full(A - w'*w*lam);
limit = 100000;

% Run the deterministic version
cBest = iterated_genlouvain(B,limit,0,0,0,cinit);
BestMod = compute_modularity(cBest,A);

% for times = 1:k
%     c = iterated_genlouvain(B,limit,0,1,0,cinit);
%      mod = compute_modularity(c,A);
%     if mod > BestMod
%         BestMod = mod;
%         cBest = c;
%     end
% end

Mods = [];
for times = 1:k
    c = iterated_genlouvain(B,limit,0,1,1,cinit);
    mod = compute_modularity(c,A);
    Mods = [Mods; mod];
    if mod > BestMod
        BestMod = mod;
        cBest = c;
    end
end



end