function A = StandardizeFully(A)
%Make graph unweighted and undirected and connected
A = spones(A|A');

[numcomps,comps] = graphconncomp(A);

B = find(comps == mode(comps));
%B = TheLargestComponent(A);
A = A(B,B);
fprintf('This graph has %d vertices and %d edges\n',size(A,2),nnz(A)/2)
A = A - diag(diag(A));
end
