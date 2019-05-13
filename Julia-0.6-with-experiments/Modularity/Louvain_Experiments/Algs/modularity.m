function Q = modularity(A,c)

Q = 0;
m = nnz(A)/2;
d = sum(A,1)';
for i = 1:max(c)
    
    
    S = find(c == i);
    for a = 1:numel(S)
        for b = a+1:numel(S)
            u = S(a);
            v = S(b);
            Q = Q + A(u,v) - d(u)*d(v)/(2*m);
        end
    end


end

n = size(A,1);
C = 0;
% Fix terms on diagonal previously not accounted for
for i = 1:n
    C = C + d(i)*d(i)/(2*m)^2;
end

Q = Q/m - C;