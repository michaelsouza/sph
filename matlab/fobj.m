function [y,g] = fobj(G, x)
g = zeros(3 * G.nnodes, 1);
y = 0;
for k = 1:G.nedges
    lij = G.l(k)^2;
    uij = G.u(k)^2;
    i = G.i(k);
    j = G.j(k);
    iid = (3*i+1):(3*i+3);
    jid = (3*j+1):(3*j+3);
    xi = x(iid);
    xj = x(jid);
    dij = norm(xi - xj)^2;

    yij = 0;
    if dij < lij
%         fprintf('%d %d %g %g\n', i, j, dij, lij);
        [yij,gij] = N(xi, xj, dij, lij);
    end
    
    if dij > uij
%         fprintf('%d %d %g %g\n', i, j, dij, uij);
        [yij,gij] = N(xi, xj, dij, uij);
    end
    
    if abs(yij) > 0
        g(iid) = g(iid) + gij;
        g(jid) = g(jid) - gij;
        y = y + yij^2;
    end
end

% fprintf('finished\n')
end

% N
function [yij, gij] = N(xi, xj, dij, bij)
    yij = (dij - bij);
    gij = 4 * yij * (xi - xj);
end