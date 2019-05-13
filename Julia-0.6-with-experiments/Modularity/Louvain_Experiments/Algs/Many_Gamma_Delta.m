function [cBest,obj] = Many_Gamma_Delta(A,Dist,gamma,delta,k)
% This is LP5. Rounding scheme to get our 5 approximation. 
%
% This repeatedly calls the gamma_delta_rounding scheme k times and returns
% the best objective value.

obj = nnz(A);
for i = 1:k
    [~,c] = gamma_delta_rounding(Dist,gamma,delta);
    
    o = CCminDisagreeObj(A,c);
    
    if o < obj
        obj = o;
        cBest = c;
    end
end