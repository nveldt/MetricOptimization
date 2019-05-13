function B = modularityLPrelax(A,D)

n = size(A,1);
d = sum(A,1)';
B = 0;
m = nnz(A)/2;
for i = 1:n
    for j = i+1:n
        B = B + (A(i,j)-d(i)*d(j)/(2*m))*(1-D(j,i));
    end
end

C = 0;
for i = 1:n
    C = C + d(i)*d(i)/(2*m)^2;
end

B = B/m - C;
end