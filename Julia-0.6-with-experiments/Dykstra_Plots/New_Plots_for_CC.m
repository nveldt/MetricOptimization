clear
Graphs = {'power', 'cagrqc', 'caHepTh', 'caHepPh'};
% Notes: 
% 1. For graph number 4, for some reason the data for the first hundred
%   iterations didn't store properly. No idea why. You can still see 
%   results for the last couple hundred iterations.
%
% 2. The constraint satisfaction plot has long flat steps because a full
%    constraint check (which takes a while) is only updated every 20
%    iterations.

for graphnum = 1:4
switch graphnum
    case 1
        n = 4941;
        ComputeTime = 7.6*3600;
    case 2
        n = 4158;
        ComputeTime = 6.6*3600;
    case 3
        n = 8638;
        ComputeTime = 88.3*3600;
    case 4
        n = 11204;
        ComputeTime = 167.5*3600;
    otherwise
        n = 1;
end


        
graphname = char(Graphs(graphnum));
load(strcat('CC_output/',char(Graphs(graphnum)),'_plot_data.mat'))
lw = 1;
gap = (duals - primals)./duals;

final = duals(end);
center = 1;

if graphnum == 4
    xs = (1:numel(duals))+93;
else
    xs = 1:numel(duals);
end

avgIter = (round(ComputeTime/numel(duals)));
ComputeTime
%% The duality gap quickly decreases, so this plot isn't as helpful

% figure(1)
% 
% % Dual objective--lower bound on OPT
% plot(xs,duals,'LineWidth',lw,'color','g')
% hold on
% 
% % Primal objective--not always an upper bound
% plot(xs,primals,'LineWidth',lw,'color','b')
% 
% % Zoom in to the important region
% top = max(max(center*primals(100:end)/final),2*final);
% ylim([-1/2*final,top])
% xlim([0, numel(duals)])
% hold off
% box off
%% It's easier to just look at the duality gap and constraint violation together

figure(2)
plot(xs,Conviolations,'LineWidth',lw,'color','k');
hold on
plot(xs,gap)
% xlim([m, numel(duals)])
xlim([xs(1),xs(end)])
legend({'Constraint Tol', 'Duality Gap'})
legend('Location','Best');
legend boxoff
%title([graphname ' (n = ' num2str(n) '), 1 iter = ' num2str(avgIter) 's']) 
box off
hold off
plotname = ['Figures/CC_' graphname]
xlabel('Number of Iterations')
set_figure_size([4,3]);
%print(plotname,'-dpdf');
print(gcf,sprintf(strcat('Figures/CC_',graphname,'.eps')),'-depsc2');
end