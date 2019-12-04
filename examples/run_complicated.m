clear; close all; clc 
addpath(genpath(fileparts(pwd)))

% make the data
nx = 2^8;
x = linspace(0,100,nx);

nt = 2^9;
t = linspace(0,20, nt);
[Tgrid,Xgrid] = meshgrid(t,x);

% wave speeds
c1 = -3*ones(size(t));
c2 = .15.*t;

u1 = exp(-.1*(Xgrid-80-c1.*Tgrid).^2);
u2 = exp(-.1*(Xgrid-c2.*Tgrid-20).^2);
u = u1+u2;

%%
params.data.u = u;
params.data.x = x.';
params.data.t = t.';

params.optim.n = 3;
params.optim.library = @(var)[ var.^2, var, ones(size(var))];
params.optim.init_lambda = 10; 

params.optim.mu = 1e8;
params.optim.zeta = 1e-4; % relaxation parameter, 1e-1
params.optim.reg = 1e4; %regularization

params.optim.type = 'ridge';

params = optimSR3(params);


params.ROM.lambda = 0.025;
params.ROM.rank = 1;

params = make_ROMS(params);

plot_complicated(params)


real = [0.15,0;0,-3;20,80];
final_err = norm(params.optim.Bsave{end}-real)/norm(real)
