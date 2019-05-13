%% Read in convergence data

fileID = fopen('CC_output/DykstraCC_output_caHepPh','r');


%% Discard the first two lines
tline = fgetl(fileID);
disp(tline)

tline = fgetl(fileID);
disp(tline)
%%

% For each line, grab the primal and dual scores and the constraint
% violation
iteration = 0;
duals = [];
primals = [];
Conviolations = [];
while ischar(tline)
    iteration = iteration + 1;
    tline = fgetl(fileID);
    
    DualPlace = strfind(tline,'Dual =');
    PrimalPlace = strfind(tline,'Primal =');
    GapPlace = strfind(tline,'gap = ');
    ConPlace = strfind(tline,'ConVio: ');
    Tri = strfind(tline,'TriVio ');
    dual = str2num(tline(DualPlace+7:PrimalPlace-3));
    primal = str2num(tline(PrimalPlace+9:GapPlace-3));
    convio = str2num(tline(ConPlace+8:Tri-3));
    
    duals = [duals; dual];
    primals = [primals; primal];
    Conviolations = [Conviolations; convio];
   
end

