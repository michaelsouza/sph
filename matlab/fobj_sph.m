function [f,g] = fobj_sph(lambda,t,G,x)
f = 0;
g = zeros(size(x));
for k = 1:G.nedges
    ik = G.i(k) + 1;
    jk = G.j(k) + 1;
    if jk > ik
        continue
    end
    lk = G.l(k);
    uk = G.u(k);
    
    iid = (3*(ik-1)+1):(3*ik);
    jid = (3*(jk-1)+1):(3*jk);
    xi = x(iid);
    xj = x(jid);
    vk = xi - xj;
    
    [f_tta_k, gk] = theta(t, vk);
    [f_ph_lk, g_ph_lk] = phi(lambda, t, lk - f_tta_k);
    [f_ph_uk, g_ph_uk] = phi(lambda, t, f_tta_k - uk);
    
    f = f + f_ph_lk + f_ph_uk;
    
%     fprintf('%d %d %g %g %g %g %g\n', ik, jk, lk, uk, f, f_ph_lk, f_ph_uk);
         
    gk = (g_ph_uk - g_ph_lk)*gk;
    g(iid) = g(iid) + gk;
    g(jid) = g(jid) - gk;
end
end