
graphs = {'power', 'cagrqc', 'caHepTh', 'caHepPh'};

% Notes: 
% 1. For graph number 4, for some reason the data for the first hundred
%   iterations didn't store properly. No idea why. You can still see 
%   results for the last couple hundred iterations.
%
% 2. The constraint satisfaction plot has long flat steps because a full
%    constraint check (which takes a while) is only updated every 20
%    iterations.

graphnum = 1;

load(strcat('CC_output/',char(graphs(graphnum)),'_plot_data.mat'))
lw = 1;
gap = (duals - primals)./duals;

final = duals(end);
center = 1;

xs = 1:numel(duals);

%% The duality gap quickly decreases, so this plot isn't as helpful

figure(1)

% Dual objective--lower bound on OPT
plot(xs,duals,'LineWidth',lw,'color','g')
hold on

% Primal objective--not always an upper bound
plot(xs,primals,'LineWidth',lw,'color','b')

% Zoom in to the important region
top = max(max(center*primals(100:end)/final),2*final);
ylim([-1/2*final,top])
xlim([0, numel(duals)])
hold off
box off
%% It's easier to just look at the duality gap and constraint violation together

figure(2)
plot(Conviolations,'LineWidth',lw,'color','k');
hold on
plot(gap)
xlim([0, numel(duals)])
legend({'Constraint Violation', 'Duality Gap'})
legend('Location','Best');
legend boxoff

box off
hold off