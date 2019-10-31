clear all, close all, clc, 
addpath(genpath(fileparts(pwd)))
sload('data/kdv.mat');

params.data.u = u;
params.data.x = x;
params.data.t = t;

params.optim.n = 2;
params.optim.library = @(var)[var.^3, var.^2, var, ones(size(var))];
params.optim.init_lambda = 2; 

params.optim.mu = 1e8;
params.optim.zeta = 1e-4; % relaxation parameter, 1e-1
params.optim.reg = 1e1; %regularization

params.optim.type = 'ridge';

params = optimSR3(params);

params.ROM.lambda = 0.01;
params.ROM.rank = 1;

params = make_ROMS(params);
plot_kdv(params)


speed = 9;
real = [0,0;0,0;speed,speed;0,-2*pi];
final_err = norm(params.optim.Bsave{end}-real)/norm(real)
