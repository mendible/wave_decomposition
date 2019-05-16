clear all, close all, clc

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
optim_params.u = u;
optim_params.x = x.';
optim_params.t = t.';

optim_params.n = 3;
optim_params.library = @(var)[ var.^2, var, ones(size(var))];
optim_params.mu = 1e8;
optim_params.zeta = 1e-4; % relaxation parameter, 1e-1
optim_params.reg = 1e4; %regularization
optim_params.lambda = 10; %0.15;
optim_params.type = 'ridge';

optim_params = optimSR3(optim_params);

lambda = 0.025;
rank = 1;

ROMS = make_ROMS(optim_params,rank,lambda);

plot_complicated(optim_params, ROMS)


% speed = 9;
% real = [0,0;0,0;speed,speed;0,0;0,-2*pi];
%
% final_err = norm(Afinal-real)/norm(real)
