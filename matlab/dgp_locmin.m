% options = optimoptions('fminunc');
function [x,f,solved] = dgp_locmin(G,ftol,display)
options = optimoptions('fminunc');
options.SpecifyObjectiveGradient = true;
options.Display = display;

x = xinit_sphere(G);
[x,f] = fminunc(@(x)fobj(G,x), x, options);
solved = is_done(G,x,ftol);
end