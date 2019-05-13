%% Runtime plot for all graphs from size 62 up to 3068


addpath('~/dev/MakeFigures/')

Times = [10; 8; 35; 21; 73; 80; 81; 166; 350; 511; 1134; 1954; 1138; 1427; 53449;25703; 34621;41080; 155333];
Ns = [62; 77; 105; 112; 115; 124; 198; 233; 297; 332; 397; 446; 453; 500;994; 1024; 1133;1222;3068];
Ms = [159; 254; 441; 425; 613; 5972; 2742; 994; 2148; 2126;914; 1413;2025;2043;3640;4916; 5451;116714;119161];

[a,b] = sort(Ms);
semilogy(Ns,Times)
loglog(Ns.*(Ns-1)/2, Times, '.-','markersize',15,'linewidth',1.25)

% figure(2)
% loglog(Ns, Times, '.-','markersize',10)

xlabel('n(n-1)/2','fontsize',10);
ylabel({'Runtime (s)'},'fontsize',10);
xlim([10^3, 10^7])
set_figure_size([2.25*1.5,1.75*1.5]);
print(gcf,sprintf('Figures/LR_runtime.eps'),'-depsc2');