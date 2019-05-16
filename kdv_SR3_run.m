clear all, close all, clc
load('kdv.mat');

u = qs;
t = t.'; 

optim_params.u = u;
optim_params.x = x;
optim_params.t = t;

optim_params.n = 2;
optim_params.library = @(var)[var.^3, var.^2, var, ones(size(var))];
optim_params.mu = 1e8;
optim_params.zeta = 1e-4; % relaxation parameter, 1e-1
optim_params.reg = 1e1; %regularization
optim_params.lambda = 2; %0.15;
optim_params.type = 'ridge';

optim_params = optimSR3(optim_params);



lambda = 0.01;
rank = 1;

ROMS = make_ROMS(optim_params,rank, lambda);
plot_kdv(optim_params,ROMS)



% opt.AB = 'B';
% opt.frames = 1:length(optim_params.Wsave);
% opt.sepFigs = 0;
% opt.fname = 'kdv.mp4';
% plot_iterations(optim_params,opt);

% speed = 9;
% real = [0,0;0,0;speed,speed;0,0;0,-2*pi];
% 
% final_err = norm(Afinal-real)/norm(real)
