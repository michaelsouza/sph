close all

%% small problem
clear
ftol = 1E-2;
display = 'off'; % final, iter

% [G,x] = read_data('/home/michael/gitrepos/bioinfo/codes/source/bp/mdgp_s0_n10_md5_dp_01.dat');
[G,x] = read_data('/home/michael/gitrepos/bioinfo/codes/source/bp/mdgp_s0_n25_md5_dp_01.dat');

[y,f,solved] = dgp_locmin(G,ftol,display);

[m,n] = size(x);
y = reshape(y,m,n);

[y,rxy] = rmsd(x,y);

fig = figure;
hold on
plot_backbone('x',x);
plot_backbone('y',y);

%% large test
clc; clear
rng(1)

ftol = 1E-2;
n = 100; d_max = 5; d_eps = 0.1; x_eps = 0.1;
display = 'off'; % final, iter
ntests = 100;

[G,x] = create_mdgp_grid_instance(n,d_max,d_eps,x_eps);
[m,n] = size(x);

fprintf(' iter\t fobj/n\t rmsd\t dme\tsolved\n');
for i = 1:ntests
    [y,f,solved] = dgp_locmin(G,ftol,display);
    y = reshape(y,m,n);
    [y,rxy] = rmsd(x,y);
    dxy = dme(x,y);
    fprintf('%5d\t%5.2g\t%5.2g\t%5.2g\t%d\n',i,f,rxy,dxy,solved);
    if solved
        y = reshape(y,m,n);
        fig = figure;
        hold on
        plot_backbone('x',x);
        plot_backbone('y',y);
    end
end