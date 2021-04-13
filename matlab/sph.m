function [x,f,solved] = sph(G,ftol)
lambda = 0.5;
rho = 0.99;
% set initial tau
tau = sort((G.l + G.u) / 2);
tau = tau(ceil(0.75 * G.nedges));

maxit = 10000;
x = xinit_sphere(G);

options = optimoptions('fminunc');
options.SpecifyObjectiveGradient = true;
options.Display = 'off';
solved = false;
fx_s = 0; gx_s = 0; flag = 0; output.funcCount = 0;
for i = 0:maxit
    f  = @(x)fobj_sph(lambda,tau,G,x);
    fx = fobj_sph(lambda,1E-16,G,x);
    if mod(i, 20) == 0
        print_header();
    end
    fprintf('%5d | %3.2e | %3.2e | %3.2e | %3.2e |   %d   | %5d\n', ...
        i, tau, fx, fx_s, norm(gx_s), flag, output.funcCount);
    
    if fx < ftol
        solved = true;
        break
    end
    [x,fx_s,flag,output,gx_s] = fminunc(f,x,options);  
    tau = rho * tau;
end
end

function print_header()
fprintf(' iter |   tau    |    fx    |  fx_s    |  |gx_s|  | flag  | #fx_s\n');
end