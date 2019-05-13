    addpath('Community_BGLL_Matlab/')
    addpath('GenLouvain-2.1/')
    addpath('Algs')
    addpath('../data')



	load Netscience
    
    [cLouRand,BestMod,Mods] = many_louvain(A,1);