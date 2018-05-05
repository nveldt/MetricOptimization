% Try all graphs with three nodes, see if the correlation clustering LP
% relaxation without constraint x_ij <= 1 ever gives a solution where some
% variable, xij, xik, or xjk is greater than one. 
%
% If not, I think that implies that we can just drop the x_ij <= 1
% constraint for all graphs



T = Get_ConstraintsFastest(3);
rhs = zeros(3,1);
%% For a full triangle, objective is xij + xik + xjk

objective = [1 1 1];

%% For a triangle with two positive edges:

objective = [1 1 -1];

%% One positive edge

objective = [1 -1 -1];

%% No positive edges

objective = [-1 -1 -1];

%% Run Gurobi
clear model
model.obj = objective; 
sense = '<';



model.A = T;
model.rhs = rhs;
model.sense = sense;
model.ub = [1;1;1];
model.vtype = 'C';
model.modelsense = 'min';

result = gurobi(model);

cond = result.objval
bset = result.x