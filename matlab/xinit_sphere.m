function x = xinit_sphere(G)
% References
% [1] Moré, J. J., & Wu, Z. (1999). Distance geometry optimization for 
%     protein structures. Journal of Global Optimization, 15(3), 219-234.
G.i = G.i + 1;
G.j = G.j + 1;
x = zeros(3 * G.nnodes, 1);
S = randperm(G.nnodes);
D = sparse(G.i,G.j,ones(G.nedges,1),G.nnodes, G.nnodes);
G.L = sparse(G.i,G.j,G.l,G.nnodes, G.nnodes);
G.U = sparse(G.i,G.j,G.u,G.nnodes, G.nnodes);
for k = 1:G.nnodes
    i = S(k);
    iid = (3*(i-1)+1):(3*(i-1)+3);
    xi  = x(iid);
    [~,j] = find(D(i,:));
    M = intersect(j,S(i+1:end));
    for j = M
        dij = (G.L(i,j) + G.U(i,j)) / 2;
        jid = (3*(i-1)+1):(3*(i-1)+3);
        % y is uniformally distributed over the S(x(j),dij)
        y = normrnd(0,1,3,1);
        y = (dij / norm(y)) * y;
        x(jid) = xi + y;
    end
end
end