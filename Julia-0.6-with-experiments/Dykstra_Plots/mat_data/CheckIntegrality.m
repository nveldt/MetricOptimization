load graphlist

% This checks to see for which graphs the convex relaxation actually
% produced the optimal scaled sparsest cut partition.
% This is the case whenever the matrix entries are either zero, or all one
% other (graph and partition dependent) constant.
for i = 1:numel(Graphs)
    
    load(strcat(char(Graphs(i)),'_output.mat'))
    [a,b,c] = find(X);
    diff = max(c) - min(c);
    
    if diff < .0001
        Graphs(i)
    end
    
end