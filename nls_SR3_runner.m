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
ROMS = make_ROMS(optim_params,rank, lambda)

colormap(flipud(gray))

figure()
surfl(x,t,sum(lowsrpca,3).'+10), shading interp, colormap(gray)
hold on 
imagesc(x,t,sum(lowsrpca,3).'/max(lowsrpca(:))), shading interp, colormap(flipud(gray))
view([12,23])
title('shifted RPCA')

figure()
surfl(x,t,sum(lowspod,3).'+10), shading interp, colormap(gray)
hold on 
imagesc(x,t,sum(lowspod,3).'/max(lowspod(:))), shading interp, colormap(flipud(gray))
view([12,23])
title('shifted POD')

figure()
surfl(x,t,sum(lowpod,3).'+10), shading interp, colormap(gray)
hold on 
imagesc(x,t,sum(lowpod,3).'/max(lowpod(:))), shading interp, colormap(flipud(gray))
view([12,23])
title('unshifted POD')


% % figure
% % subplot(1,3,1)
% % pcolor(u.'), shading interp, colormap(flipud(gray)), axis off
% % title('Original Data','fontsize',18)
% % subplot(1,3,2)
% % pcolor(unew(:,:,1).'), shading interp, colormap(flipud(gray)), axis off
% % title('Shifted for Right Wave','fontsize',18)
% % subplot(1,3,3)
% % pcolor(unew(:,:,2).'), shading interp, colormap(flipud(gray)), axis off
% % title('Shifted for Left Wave','fontsize',18)
% % set(gcf,'Position',[439, 1475, 1421, 357])

% figure
% subplot(2,2,1)
% imagesc(ulow(:,:,1))
% subplot(2,2,3)
% imagesc(ulow(:,:,2))
% subplot(2,2,2)
% imagesc(unew2(:,:,1))
% subplot(2,2,4)
% imagesc(unew2(:,:,2))

% figure
% subplot(1,3,2)
% pcolor(u(:,:,1).'), shading interp, colormap(flipud(gray)), axis off
% title('Right Wave, 10 mode reconstruction','fontsize',18)
% subplot(1,3,3)
% pcolor(lowsrpca(:,:,2).'), shading interp, colormap(flipud(gray)), axis off
% title('Left Wave, 10 mode reconstruction','fontsize',18)
% subplot(1,3,1)
% pcolor(sum(lowsrpca,3).'), shading interp, colormap(flipud(gray)), axis off
% title('20 mode reconstruction','fontsize',18)
% set(gcf,'Position',[439, 1475, 1421, 357])



%%
%
% data = figure(2);
% w = waterfall(t(1:5:end),x,u(:,1:5:end));
% set(w, 'EdgeColor', [0 0 0]);
% set(gca,'fontsize', 20, 'Xtick', [-3, 3], 'Ytick', [0, 0.4], 'Ztick', [40]);
% view([20, 75]);
% xlim([-pi, pi]);

% opt.AB = 'B';
% opt.sepFigs = 0;
% opt.frames = 1:100:length(optim_params.Wsave);
% % opt.fname = 'nls.mp4';
% plot_iterations(optim_params,opt);

