clear all, close all, clc
load('burgers.mat');

u = abs(u);
x = x.';
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
shifts = optim_params.shifts;
%%
dx = x(2)-x(1);
[n1,n2] = size(u);
mu = n1*n2/(4*sum(abs(u(:))));
lambda = 1/sqrt(max(n1,n2));
rank = 1;

[U.pod,S.pod,V.pod] = svd(u);
lowpod = U.pod(:,1:rank)*S.pod(1:rank,1:rank)*V.pod(:,1:rank)';

for j = 1:size(shifts,2)
    %     [~,tmp] = min(abs(abs(shifts(:,j))-repmat(t',[1 nt])'));
    %     shift_ind(:,j) = tmp.*sign(shifts(:,j))';
    for jj = 1:length(t)
        unew(:,jj,j) = circshift(u(:,jj),-round(shifts(jj,j)/dx));
    end
    %     [LowRank{j},Sparse{j}] = RPCA(unew(:,:,j),lambda,mu);
    [U.spod{j},S.spod{j},V.spod{j}] = svd(unew(:,:,j));
    [LowRank{j},Sparse{j},rpcaiter] = inexact_alm_rpca(unew(:,:,j),0.01);
    [U.srpca{j},S.srpca{j},V.srpca{j}] = svd(LowRank{j});
   
    usrpca(:,:,j) = U.srpca{j}(:,1:rank)*S.srpca{j}(1:rank,1:rank)*V.srpca{j}(:,1:rank)';
    uspod(:,:,j) = U.spod{j}(:,1:2*rank)*S.spod{j}(1:2*rank,1:2*rank)*V.spod{j}(:,1:2*rank)';
    
    for jj = 1:length(t)
        lowsrpca(:,jj,j) = circshift(usrpca(:,jj,j),round(shifts(jj,j)/dx));
        lowspod(:,jj,j) = circshift(uspod(:,jj,j),round(shifts(jj,j)/dx));
    end
end

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

%%

% % data = figure(2);
% % w = waterfall(x,t(1:5:end),u(:,1:5:end).');
% % set(w, 'EdgeColor', [0 0 0]);
% % set(gca,'fontsize', 20, 'Xtick', [-3, 3], 'Ytick', [0, 0.4], 'Ztick', [40]);
% % view([20, 75]);
% % xlim([-pi, pi]);
% 
% opt.AB = 'B';
% opt.frames = 1:length(optim_params.Wsave);
% opt.sepFigs = 0;
% opt.fname = 'burgers.mp4';
% plot_iterations(optim_params,opt);

