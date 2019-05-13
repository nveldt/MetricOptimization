% tag = 'gamma' means the results with gamma = 5.0
% tag = '20.0' means small experiments where we used 20
% tag = 2.0 means large experiments where we used 2

% This script file in particular shows modularity bounds and times for
% five small graphs whose exact modularity is known. The bounds are very
% good, and the times are competitive for 


addpath('Community_BGLL_Matlab/')
addpath('GenLouvain-2.1/')
addpath('Algs')
addpath('../output/')
addpath('~/GitHubRepos/Tri-Con-LP/Code/graphs/')

gammatag = '20.0';
graphs = {'KarateA','dolphinsA','lesmisA','polbooksA','footballA'};
 
 for gr = 1:5
    graph = graphs(gr);
    
    load(char(graph));

    %% Load the data
    load(char(strcat(graph,'_Gam_',char(gammatag),'_DykstraModularity_output.mat')))

    delta = Aposteriori-1;
    D = tril(D);
    D = D+D';

    %% Get an upper bound on modularity from the D matrix
    


    % Evaluate the modularity LP relaxaction
    Bound = modularityLPrelax(A,D);
    g = 0;
    d = sum(A,2);
    m = nnz(A)/2;
    [I,J] = find(triu(A));
    for t = 1:numel(I)
        g = g+ d(I(t))*d(J(t))/(2*m);
    end
    n = size(A,1);
    C = 0;
    % Fix terms on diagonal previously not accounted for
    for i = 1:n
        C = C + d(i)*d(i)/(2*m)^2;
    end
    
    Breal = Bound/(1+delta) + delta/(1+delta)*(1 - g/m + C);

    fprintf('Modularity bound: %s, %f, %f\n',char(graph),Breal,DykstraTime)
    
 end