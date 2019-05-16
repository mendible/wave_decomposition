clear all, close all, clc
load('burgers.mat');

u = abs(u);

optim_params.u = u;
optim_params.x = x;
optim_params.t = t;

optim_params.n = 1;
optim_params.library = @(var)[sqrt(abs(var)), var, ones(size(var))];
optim_params.mu = 1e8;
optim_params.zeta = 1e-4; % relaxation parameter, 1e-1
optim_params.reg = 1e1; %regularization
optim_params.lambda = 0.1; %0.15;
optim_params.type = 'edge';

optim_params = optimSR3(optim_params);

lambda = 0.01;
rank = 1;

ROMS = make_ROMS(optim_params,rank, lambda);

plot_burgers(optim_params, ROMS)



