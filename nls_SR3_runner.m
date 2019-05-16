clear all, close all, clc
load('nls.mat');

u = abs(usol).';

x = x.';
optim_params.u = u;
optim_params.x = x;
optim_params.t = t;

optim_params.n = 3;
optim_params.library = @(var)[var.^3, var.^2, var, ones(size(var))];
optim_params.mu = 1e8;
optim_params.zeta = 1e-4; % relaxation parameter, 1e-1
optim_params.reg = 1e1; %regularization
optim_params.lambda = 0.0001; %0.15;
optim_params.type = 'ridge';

optim_params = optimSR3(optim_params);

lambda = 0.01; 
rank = 1;
ROMS = make_ROMS(optim_params,rank, lambda);

plot_nls(optim_params, ROMS)
