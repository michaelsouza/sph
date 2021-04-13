function ok = is_done(G,x,ftol)
for k = 1:G.nedges
    i = G.i(k) + 1;
    j = G.j(k) + 1;
    lij = G.l(k);
    uij = G.u(k);
    iid = (3*(i-1)+1):(3*(i-1)+3);
    jid = (3*(j-1)+1):(3*(j-1)+3);
    xi = x(iid);
    xj = x(jid);
    dij = norm(xi - xj);
    ok = lij * (1 - ftol) < dij && dij < uij * (1 + ftol);
    if ~ok 
        return
    end
end
end