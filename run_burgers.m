clear all, close all, clc
load('data/burgers.mat');

params.data.u = u;
params.data.x = x;
params.data.t = t;

params.optim.n = 1;
params.optim.library = @(var)[sqrt(abs(var)), var, ones(size(var))];
params.optim.init_lambda = 0.1; 

params.optim.zeta = 1e-4; % relaxation parameter, 1e-1
params.optim.reg = 1e1; %regularization

params.optim.type = 'edge';

params = optimSR3(params);

params.ROM.lambda = 0.01;
params.ROM.rank = 1;

params = make_ROMS(params);

plot_burgers(params)



