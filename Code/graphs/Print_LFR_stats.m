%% Print some stats for LFR graphs

for n = [200 300 500 1000 1500 2000 2500]
    
    load(strcat('A',num2str(n)))
    for lam = [.9]
     fprintf('LFR%d & %d & %d & %.2f &  \\\\ \n',n,n,nnz(A)/2, lam)    
    end
    for lam = [.75 .5 .25 .1 .05 .01]
    fprintf(' &  &  & %.2f &  \\\\ \n', lam)
    end
    fprintf('\\midrule\n')
end