clear all, close all, clc

nx = 2^8;
x = linspace(0,100,nx);

nt = 2^9;
t = linspace(0,20, nt);
[Tgrid,Xgrid] = meshgrid(t,x);

c1 = -3*ones(size(t));
c2 = .15.*t;

u1 = exp(-.1*(Xgrid-80-c1.*Tgrid).^2);
u2 = exp(-.1*(Xgrid-c2.*Tgrid-20).^2);
u = u1+u2;

%%

% filter1 = exp(-.05*(Xgrid-80).^2);
% filter2 = exp(-.05*(Xgrid-20).^2);
% filter1 = filter1>0.05;
% filter2 = filter2>0.05;
% plot(u(:,1),'k')
% hold on
% plot(filter1(:,1),'r')
% plot(filter2(:,1),'b')
% plot(u(:,1).*filter1(:,1),'y')
% plot(u(:,1).*filter2(:,1),'g')

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

    for jj = 1:length(t)
        unew(:,jj,j) = circshift(u(:,jj),-round(shifts(jj,j)/dx));
    end
    [U.spod{j},S.spod{j},V.spod{j}] = svd(unew(:,:,j));
    [LowRank{j},Sparse{j},rpcaiter] = inexact_alm_rpca(unew(:,:,j),0.025);
    [U.srpca{j},S.srpca{j},V.srpca{j}] = svd(LowRank{j});
   
    usrpca(:,:,j) = U.srpca{j}(:,1:rank)*S.srpca{j}(1:rank,1:rank)*V.srpca{j}(:,1:rank)';
    uspod(:,:,j) = U.spod{j}(:,1:2*rank)*S.spod{j}(1:2*rank,1:2*rank)*V.spod{j}(:,1:2*rank)';
    
    for jj = 1:length(t)
        lowsrpca(:,jj,j) = circshift(usrpca(:,jj,j),round(shifts(jj,j)/dx));
        lowspod(:,jj,j) = circshift(uspod(:,jj,j),round(shifts(jj,j)/dx));
    end
end


%%
figure(1)
surfl(x,t,u.'+10), shading interp, colormap(gray)
hold on 
imagesc(x,t,u.'/max(u(:))), shading interp, colormap(flipud(gray))
view([12,23])
title('Original Data')

figure(2)
surfl(x,t,(sum(lowsrpca,3)/2).'+10), shading interp, colormap(gray)
hold on 
imagesc(x,t,sum(lowsrpca,3).'/max(lowsrpca(:))), shading interp, colormap(flipud(gray))
view([12,23])
title('Shifted RPCA')

figure(3)
surfl(x,t,(sum(lowspod,3)/2).'+10), shading interp, colormap(gray)
hold on 
imagesc(x,t,sum(lowspod,3).'/max(lowspod(:))), shading interp, colormap(flipud(gray))
view([12,23])
title('Shifted POD')

figure(4)
surfl(x,t,(sum(lowpod,3)/2).'+10), shading interp, colormap(gray)
hold on 
imagesc(x,t,sum(lowpod,3).'/max(lowpod(:))), shading interp, colormap(flipud(gray))
view([12,23])
title('Unshifted POD')


%%

% data = figure();
% w = waterfall(x,t(1:5:end),u(:,1:5:end).');
% set(w, 'EdgeColor', [0 0 0]);
% set(gca,'fontsize', 20, 'Xtick', [-3, 3], 'Ytick', [0, 0.4], 'Ztick', [40]);
% view([20, 75]);
% xlim([-pi, pi]);
% % saveas(data,'original_kdv.png')
% 
% 
% opt.AB = 'B';
% opt.frames = [1:100:length(optim_params.Wsave)];
% % opt.frames = 1;
% % opt.frames = length(optim_params.Wsave);
% opt.sepFigs = 0;
% % opt.fname = 'complicated_reg1e4_z1e-1.mp4';
% plot_iterations(optim_params,opt);

% speed = 9;
% real = [0,0;0,0;speed,speed;0,0;0,-2*pi];
%
% final_err = norm(Afinal-real)/norm(real)
