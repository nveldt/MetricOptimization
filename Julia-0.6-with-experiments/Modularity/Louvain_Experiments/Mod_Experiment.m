function Mod_Experiment(graph,gammatag)

    % tag = 'gamma' means the results with gamma = 5.0
    % tag = '20.0' means small experiments where we used 20
    % tag = 2.0 means large experiments where we used 2
    addpath('Community_BGLL_Matlab/')
    addpath('GenLouvain-2.1/')
    addpath('Algs')
    addpath('../output/')
    addpath('~/GitHubRepos/Tri-Con-LP/Code/graphs/')
    fid = fopen(strcat(graph,gammatag,'_Louvain_results.txt'),'w');

    load(graph);

    %% Load the data
    load(strcat(graph,'_Gam_',gammatag,'_DykstraModularity_output'))

    delta = Aposteriori-1;
    D = tril(D);
    D = D+D';

    %% Get an upper bound on modularity from the D matrix
    
    n = size(A,1);
    C = 0;
    % Fix terms on diagonal previously not accounted for
    for i = 1:n
        C = C + d(i)*d(i)/(2*m)^2;
    end

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

    fprintf(fid,'Modularity bound: %f \n',Breal)

    lTimes = 15;

    %% Randomized louvain for better results

    [cLouRand,BestMod,Mods] = many_louvain(A,lTimes);

    fprintf(fid,'Randomized Louvain: Best = %f and Median = %f \n',BestMod, median(Mods))

    %% CGW rounding, with refinement

    gamma = 1/2;
    delta = 1/4;
    Basics = [];
    refs = [];
    for times = 1:lTimes
        k = 50;
        [cCGW,~] = Many_Gamma_Delta(A,D,gamma,delta,k);
        Qcgw = modularity(A,cCGW);

        Basics = [Basics; Qcgw];
        [cRef, Qcgw_r] = many_louvain(A,1,cCGW);
        refs = [refs; Qcgw_r];
    end

    fprintf(fid,'CGW: %f \t %f \n',max(Basics), median(Basics))
    fprintf(fid,'CGW + Louv: %f \t %f \n',max(refs), median(refs))


    %% ThreeLP rounding, cutoff = 1/2

    Basics = [];
    refs = [];

    for time = 1:lTimes
        k = 50;
        [c3LP,obj]= ThreeLP_round(A,D,k,1/2);
        Q3lp = modularity(A,c3LP);
        [cRef,Q3lp_r] = many_louvain(A,1,c3LP);
        Basics = [Basics; Q3lp];
        refs = [refs;Q3lp_r];

    end

    fprintf(fid,'3LP 1: %f \t %f \n',max(Basics), median(Basics))
    fprintf(fid,'3LP 1 + Louv: %f \t %f \n',max(refs), median(refs))

    %% ThreeLP rounding, cutoff = 1/3

    Basics = [];
    refs = [];

    for time = 1:lTimes
        k = 50;
        [c3LP,obj]= ThreeLP_round(A,D,k,1/3);
        Q3lp = modularity(A,c3LP);
        [cRef,Q3lp_r] = many_louvain(A,1,c3LP);
        Basics = [Basics; Q3lp];
        refs = [refs;Q3lp_r];

    end

    fprintf(fid,'3LP 2: %f \t %f \n',max(Basics), median(Basics))
    fprintf(fid,'3LP 2 + Louv: %f \t %f \n',max(refs), median(refs))

end
