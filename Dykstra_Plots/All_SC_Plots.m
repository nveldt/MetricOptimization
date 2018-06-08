%% Convergence plots for Running DykstraSC

clear
addpath('mat_data')
load graphlist

% Numbering each graph for easy visualization:
% 
% 1. Jazz
% 2. SmallW
% 3. C. Elegans Neural
% 4. USAir97
% 5. Netscience
% 6. Erdos991
% 7. C. Elegans Metabolic
% 8. Harvard500
% 9. Roget
% 10. SmaGri
% 11. Email
% 12. Polblogs
% 13. Vassar85
% 
% Smaller Graphs not included in the paper:
% 
% 14. Les Mis
% 15. Dolphins
% 16. Polbooks
% 17. Football
% 18. Adj-noun

%% select a graph

for graphnum = 1:13
load(strcat('mat_data/',char(Graphs(graphnum)),'_output.mat'))
graphname = char(Graphs(graphnum));
if strcmp(graphname(end),'A')
    graphname = graphname(1:end-1);
end
%graphnum
n = size(X,1);
avgIter = round(ComputeTime/numel(duals),2);

fprintf(' case %d \n \t title( %s (n = %d), 1 iter \\approx %.1f s) \n',graphnum,graphname,n,avgIter)


%% Scale the Primal and Dual Objectives
% Then you can visualize the Primal/Dual gap and constraint satisfaction at
% the same time.


lw = 1;
final = duals(end);
xs = 1:numel(duals);
center = .1;

plotname = ['Figures/' graphname];

figure(1)

% Scaled/normalized dual objective--lower bound on normalized OPT
plot(xs,duals/final*center,'LineWidth',lw,'color','g')
hold on

% Scaled/normalized primal objective--not always an upper bound
plot(xs,primals/final*center,'LineWidth',lw,'color','b')

% Dotted line indicating optimal value
plot(xs,center*ones(numel(primals),1),'--','LineWidth',lw,'color','r')

% Maximum constraint violation
plot(Conviolation,'LineWidth',lw,'color','k');

% Zoom in to the important region
top = max(max(center*primals(100:end)/final),2*center);
ylim([-1/2*center,top])
xlim([50, numel(duals)])

legend({'0.1*Dual/OPT', '0.1*Primal/OPT','0.1 = Scaled OPT','Constraint Tol'})
legend('Location','Best');
legend boxoff
xlabel('Number of Iterations')
ylabel('Constraint Tol/ QP Scores')
switch graphnum
 case 1 
 	 title( 'Jazz (n = 198), 1 iter = 0.1 s') 
 case 2 
 	 title( 'SmallW (n = 233), 1 iter = 0.1 s') 
 case 3 
 	 title( 'C.El. Neural (n = 297), 1 iter = 0.2 s') 
 case 4 
 	 title( 'USair97 (n = 332), 1 iter = 0.3 s') 
 case 5 
 	 title( 'Netscience (n = 379), 1 iter = 0.2 s') 
 case 6 
 	 title( 'Erdos991 (n = 446), 1 iter = 0.3 s') 
 case 7 
 	 title( 'C.El. Metabolic (n = 453), 1 iter = 0.4 s') 
 case 8 
 	 title( 'Harvard500 (n = 500), 1 iter = 0.6 s') 
 case 9 
 	 title('Roget (n = 994), 1 iter = 1.8 s') 
 case 10 
 	 title( 'SmaGri (n = 1024), 1 iter = 1.5 s') 
 case 11 
 	 title( 'Email (n = 1133), 1 iter = 2.0 s') 
 case 12 
 	 title( 'Polblogs (n = 1222), 1 iter = 4.9 s') 
 case 13 
 	 title( 'Vassar85 (n = 3068), 1 iter = 52 s') 
    otherwise
        title([graphname ' (n = ' num2str(n) '), 1 iter = ' num2str(avgIter) 's'])    
end
box off
hold off

set_figure_size([4,3]);
print(plotname,'-dpdf');

end