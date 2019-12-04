clear; close all; clc 
addpath(genpath(fileparts(pwd)))
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

params.ROM.lambda = 0.1;
params.ROM.rank = 1;

params = make_ROMS(params);

plot_burgers(params)


xind = params.optim.xind;
tind = params.optim.tind;
for jj = 1:length(xind)
    uamp(jj) = u(xind(jj),tind(jj));
end
tpts = params.optim.tpts;

deriv = params.optim.Bsave{end}(1)./(2*sqrt(abs(tpts)));

figure
yy = smooth(uamp,50);
plot(deriv, 'linewidth', 3), hold on, plot(yy, 'linewidth', 3)
xticks([])
yticks([])


