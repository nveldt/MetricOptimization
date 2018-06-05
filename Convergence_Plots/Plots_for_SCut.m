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

graphnum = 8;
load(strcat('mat_data/',char(Graphs(graphnum)),'_output.mat'))

%% Plot the primal and dual objective scores

n = size(X,1);
lw = 1;
final = duals(end);
xs = 1:numel(duals);

figure(1)
plot(xs,duals,'LineWidth',lw,'color','g')
hold on
plot(xs,primals,'LineWidth',lw,'color','b')
plot(xs,final*ones(numel(primals),1),'--','LineWidth',lw,'color','r')

% Zoom in to the important region
top = max(max(primals(100:end)),2*final);
ylim([-1/2*final,top])
xlim([50, numel(duals)])

%% Plot the maximum constraint violation

figure(2)
plot(xs,Conviolation)
xlim([50, numel(duals)])
ylim([0,1])

%% The duality gap isn't that helpful to look at
% It only imparts useful information when k is large enough so that the
% primal objective is an upper bound on OPT, which often it isn't.

figure(3)
plot(xs, (duals-primals)./duals)

%% Scale the Primal and Dual Objectives
% Then you can visualize the Primal/Dual gap and constraint satisfaction at
% the same time.

center = .1;


figure(graphnum+3)

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

legend({'0.1*Dual/OPT', '0.1*Primal/OPT','0.1 = 0.1*OPT/OPT','Constraint Tol'})
legend('Location','Best');
legend boxoff
title(strcat(char(Graphs(graphnum)),'\_',num2str(n)))

box off
hold off