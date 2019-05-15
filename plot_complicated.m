function plot_complicated(params, ROMS)
close all

Asave = params.Asave;
Bsave = params.Bsave;
Wsave = params.Wsave;
xpts = params.xpts;
tpts = params.tpts;
t = params.t;
x = params.x;
u = params.u;
n = size(Asave{1},2);

colors_mat = {[0,0.4470,0.7410],...
    [0.8500,0.3250,0.0980],...
    [0.9290,0.6940,0.1250],'r','k','c'};

L = params.library(t);

usrpca = ROMS.usrpca;
uspod = ROMS.uspod;
upod = ROMS.upod;

xlims = [-0.2232 99.7768];
tlims = [0.0595 20.0595];
views = [12,23];

xpos_flat = [-5 10 -1];
ypos_flat = [50 -0.5 0];

xpos = [50 -2 0];
ypos = [125 10 0];
zpos = [-27 0 1.5000];

f1 = figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8]);
% plot the 3D surface of the data
surfl(x,t,u.'+2), shading interp, colormap(gray)
hold on
pcolor(x,t,u.'), shading interp, colormap(flipud(gray))
% title('Original Data','fontsize',18)
set_3d_figs(views, xpos, ypos, zpos)


f2 = figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8]);
pcolor(x,t,u.'), shading interp, colormap(flipud(gray))
% title('Original Data','fontsize',18)
set_flat_figs(xlims,tlims,xpos_flat,ypos_flat)

f3 = figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8]);
% plot the edge/ridge-found points
pcolor(x,t,u.'), shading interp, colormap(flipud(gray))
hold on
scatter(xpts,tpts,'c','filled')
% title('Ridge Detection','fontsize',18)
set_flat_figs(xlims,tlims,xpos_flat,ypos_flat)


% plot initial spectral clustering
f4 = figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8]);
pct = sum(Wsave{1})/sum(Wsave{1}(:));
pcolor(x,t,u.'), shading interp, colormap(flipud(gray))
hold on
for jj = 1:n
    if pct(jj) > 0.01
        for kk = 1:N
            if Wsave{1}(kk,jj) == 1
                plot(xpts(kk),tpts(kk),'.','markersize',16,'color',colors_mat{jj})
            end
        end
    end
end
% title('Spectral Clustering','fontsize',18)
set_flat_figs(xlims,tlims,xpos_flat,ypos_flat)


f5 = figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8]);
% plot initial models
shifts = L*Bsave{1};
pcolor(x,t,u.'), shading interp, colormap(flipud(gray))
hold on
for jj = 1:n
    if pct(jj) > 0.01
        xtmp = shifts(:,jj);
        plot(xtmp,t,'color',colors_mat{jj},'LineWidth',6)
    end
end
% title('Model Discovery','fontsize',18)
set_flat_figs(xlims,tlims,xpos_flat,ypos_flat)

% plot final spectral clustering
f6 = figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8]);
pct = sum(Wsave{end})/sum(Wsave{end}(:));
pcolor(x,t,u.'), shading interp, colormap(flipud(gray))
hold on
for jj = 1:n
    if pct(jj) > 0.01
        for kk = 1:N
            if Wsave{end}(kk,jj) == 1
                plot(xpts(kk),tpts(kk),'.','markersize',16,'color',colors_mat{jj})
            end
        end
    end
end
% title(['Final Clusters, Iteraton ',num2str(length(Asave))],'fontsize',18)
set_flat_figs(xlims,tlims,xpos_flat,ypos_flat)


% plot final models
f7 = figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8]);
shifts = L*Bsave{end};
pcolor(x,t,u.'), shading interp, colormap(flipud(gray))
hold on
for jj = 1:n
    if pct(jj) > 0.01
        xtmp = shifts(:,jj);
        plot(xtmp,t,'color',colors_mat{jj},'LineWidth',6)
    end
end
% title(['Final Models, Iteraton ',num2str(length(Asave))],'fontsize',18)
set_flat_figs(xlims,tlims,xpos_flat,ypos_flat)


%%

f8 = figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8]);
surfl(x,t,usrpca.'+2), shading interp, colormap(gray)
hold on 
imagesc(x,t,usrpca.'/max(usrpca(:))), shading interp, colormap(flipud(gray))
set_3d_figs(views, xpos, ypos, zpos)
title('Shifted RPCA')

f9 = figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8]);
surfl(x,t,uspod.'+2), shading interp, colormap(gray)
hold on 
imagesc(x,t,uspod.'/max(uspod(:))), shading interp, colormap(flipud(gray))
set_3d_figs(views, xpos, ypos, zpos)
title('Shifted POD')

f10 = figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8]);
surfl(x,t,upod.'+2), shading interp, colormap(gray)
hold on 
imagesc(x,t,upod.'/max(upod(:))), shading interp, colormap(flipud(gray))
set_3d_figs(views, xpos, ypos, zpos)
title('Unshifted POD')

% 
% print(f1,'figures/comp_data','-depsc2')
% print(f2,'figures/comp_data_flat','-djpeg','-r600')
% print(f3,'figures/comp_init_ridge','-depsc2')
% print(f4,'figures/comp_init_clusters','-depsc2')
% print(f5,'figures/comp_init_models','-depsc2')
% print(f6,'figures/comp_clusters','-depsc2')
% print(f7,'figures/comp_models','-depsc2')
% print(f8,'figures/comp_srpcs','-depsc2')
% print(f9,'figures/comp_spod','-depsc2')
% print(f10,'figures/comp_pod','-depsc2')


end

