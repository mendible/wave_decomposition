clear all, close all, clc
load('data/kdv.mat');

params.data.u = u;
params.data.x = x;
params.data.t = t;

params.optim.n = 2;
params.optim.library = @(var)[var.^3, var.^2, var, ones(size(var))];
params.optim.mu = 1e8;
params.optim.zeta = 1e-4; % relaxation parameter, 1e-1
params.optim.reg = 1e1; %regularization
params.optim.lambda = 2; %0.15;
params.optim.type = 'ridge';

params = optimSR3(params);

params.ROM.lambda = 0.01;
params.ROM.rank = 1;

params = make_ROMS(params);
plot_kdv(params)



% opt.AB = 'B';
% opt.frames = 1:length(optim_params.Wsave);
% opt.sepFigs = 0;
% opt.fname = 'kdv.mp4';
% plot_iterations(optim_params,opt);

% speed = 9;
% real = [0,0;0,0;speed,speed;0,0;0,-2*pi];
% 
% final_err = norm(Afinal-real)/norm(real)
