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
addpath('../graphs/')
addpath('~/GitHubRepos/Tri-Con-LP/Code/graphs/')

% gammatag = '20.0';
% graphs = {'KarateA','dolphinsA','lesmisA','polbooksA','footballA'};

gammatag = '2.0';
graphs = {'Netscience','SmaGriA','emailA','polblogsA','Vassar85'};
graphs = {'Netscience','polblogsA','emailA'}
gammatag = 'gamma'
gammatag = '2.0'
graphs = {'RogetA','Erdos991A','celegansmetabolicA','Harvard500A'}


 for gr = 1:numel(graphs)
    graph = graphs(gr);
    
    load(char(graph));

    %% Load the data
    load(char(strcat(graph,'_Gam_',char(gammatag),'_DykstraModularity_output.mat')))

    % Output results to a file
    fid = fopen(char(strcat('TableResults/',graph,'_',gammatag,'_Louvain_results.txt')),'w');

    delta = Aposteriori-1;
    D = tril(D);
    D = D+D';

    %% Get an upper bound on modularity from the D matrix
    
    % Evaluate the modularity LP relaxation
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
    
    % Bound from paper
    Mbound = Bound/(1+delta) + delta/(1+delta)*(1 - g/m + C);

    fprintf(fid,'Modularity bound: %s, %f, %f\n',char(graph),Mbound,DykstraTime);
    
    lTimes = 15;

    %% Randomized louvain for better results

    [cLouRand,BestModLouv,ModsLouv] = many_louvain(A,lTimes);

    fprintf(fid,'Randomized Louvain: Best = %f and Median = %f \n',BestModLouv, median(ModsLouv));

    %% CGW rounding, with refinement

    gamma = 1/2;
    delta = 1/4;
    BasicsCGW = [];
    refCGW = [];
    for times = 1:lTimes
        [cCGW,~] = Many_Gamma_Delta(A,D,1/2,1/4,50);
        Qcgw = modularity(A,cCGW);

        BasicsCGW = [BasicsCGW; Qcgw];
        [cRef, Qcgw_r] = many_louvain(A,1,cCGW);
        refCGW = [refCGW; Qcgw_r];
    end

    fprintf(fid,'CGW: %f \t %f \n',max(BasicsCGW), median(BasicsCGW));
    fprintf(fid,'CGW + Louv: %f \t %f \n',max(refCGW), median(refCGW));



    %% ThreeLP rounding, cutoff = 1/3
    Basics3LP = [];
    ref3LP = [];

    for time = 1:lTimes
        [c3LP,obj]= ThreeLP_round(A,D,50,1/3);
        Q3lp = modularity(A,c3LP);
        [cRef,Q3lp_r] = many_louvain(A,1,c3LP);
        Basics3LP = [Basics3LP; Q3lp];
        ref3LP = [ref3LP;Q3lp_r];

    end

    fprintf(fid,'3LP 2: %f \t %f \n',max(Basics3LP), median(Basics3LP));
    fprintf(fid,'3LP 2 + Louv: %f \t %f \n',max(ref3LP), median(ref3LP));

    fprintf(fid,'%s & %d & %.4f & \\emph{best}: & %.4f & %.4f & %.4f & %.4f & %.4f\n',char(graph),round(DykstraTime),Mbound,BestModLouv,max(BasicsCGW),max(refCGW),max(Basics3LP),max(ref3LP));
    fprintf(fid,'$n = %d $ &  &  & \t \\emph{med}: & %.4f & %.4f & %.4f & %.4f & %.4f\n',size(A,1),median(ModsLouv),median(BasicsCGW),median(refCGW),median(Basics3LP),median(ref3LP));

    fprintf(fid,'Extended Results\n');
    fprintf(fid,'%s & %d & %.4f & \\emph{best}: & %.4f & %.4f & %.4f \n',char(graph),round(DykstraTime),Mbound,BestModLouv,max(Basics3LP),max(ref3LP));
    fprintf(fid,'$n = %d $ &  &  & \t \\emph{med}: & %.4f & %.4f & %.4f  \n',size(A,1),median(ModsLouv),median(Basics3LP),median(ref3LP));

        fprintf('%s & %d & %.4f & \\emph{best}: & %.4f & %.4f & %.4f & %.4f & %.4f\n',char(graph),round(DykstraTime),Mbound,BestModLouv,max(BasicsCGW),max(refCGW),max(Basics3LP),max(ref3LP));
    fprintf('$n = %d $ &  &  & \t \\emph{med}: & %.4f & %.4f & %.4f & %.4f & %.4f\n',size(A,1),median(ModsLouv),median(BasicsCGW),median(refCGW),median(Basics3LP),median(ref3LP));

    fprintf('Extended Results\n');
    fprintf('%s & %d & %.4f & \\emph{best}: & %.4f & %.4f & %.4f \n',char(graph),round(DykstraTime),Mbound,BestModLouv,max(Basics3LP),max(ref3LP));
    fprintf('$n = %d $ &  &  & \t \\emph{med}: & %.4f & %.4f & %.4f  \n',size(A,1),median(ModsLouv),median(Basics3LP),median(ref3LP));

    
    fclose(fid);
    
 end