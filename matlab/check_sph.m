close all

%% small problem
clear; clc
ftol = 1E-2;

rng(1)

% [G,x] = read_data('/home/michael/gitrepos/bioinfo/codes/source/bp/mdgp_s0_n10_md5_dp_01.dat');
% [G,x] = read_data('/home/michael/gitrepos/bioinfo/codes/source/bp/mdgp_s0_n10_md5_dp_01.dat');

[y,f,solved] = sph(G,ftol);
if ~solved
    error('Small instance could not be solved')
end

[m,n] = size(x);
y = reshape(y,m,n);

[y,rxy] = rmsd(x,y);

fig = figure;
hold on
plot_backbone('x',x);
plot_backbone('y',y);
title(sprintf('rmsd = %f',rxy));

%% large test
clc; clear
rng(1)

ftol = 1E-2;
n = 100; d_max = 5; d_eps = 0.1; x_eps = 0.1;
ntests = 100;

[G,x] = create_mdgp_grid_instance(n,d_max,d_eps,x_eps);
[m,n] = size(x);

solved = zeros(ntests,1);
rxy = solved;
dxy = solved;
fx  = solved;
for i = 1:ntests
    [y,fx,solved(i)] = sph(G,ftol);
    y = reshape(y,m,n);
    [y,rxy(i)] = rmsd(x,y);
    dxy(i) = dme(x,y);
end

fprintf(' iter\t fobj/n\t rmsd\t dme\tsolved\n');
for i = 1:ntests
    fprintf('%5d\t%5.2e\t%5.2e\t%5.2e\t%d\n',i,fx(i),rxy(i),dxy(i),solved(i));
end