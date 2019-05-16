function plot_nls(params, ROMS)
close all

Asave = params.SR3.Asave;
Bsave = params.SR3.Bsave;
Wsave = params.SR3.Wsave;
xpts = params.SR3.xpts;
tpts = params.SR3.tpts;

t = params.data.t;
x = params.data.x;
u = params.data.u;
N = params.data.N;
n = params.optim.n;

colors_mat = {[1, 123, 118]/255,...
    [0.9290,0.6940,0.1250],...
    [255, 82, 0]/255};

L = params.optim.library(t);
pct = sum(Wsave{1})/sum(Wsave{1}(:));

usrpca = params.ROM.usrpca;
uspod = params.ROM.uspod;
upod = params.ROM.upod;

views = [12,23];
xpos =  [0 -0.8 0];
ypos =  [17 pi 0];
zpos = [-18 0 3.5];
zlims = [0 7.3];

xlims = [-15.2679 14.7321];
ylims = [0.0187 6.3019];
xpos_flat = [-17 3.1603 0];
ypos_flat = [-0.2678 -0.3085 0];

upshift = 3;

f1 = figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8]);
% plot the 3D surface of the data
surfl(x,t,u.'+upshift), shading interp, colormap(gray)
hold on
pcolor(x,t,u.'/max(u(:))), shading interp, colormap(flipud(gray))
% title('Original Data','fontsize',18)
set_3d_figs(views, xpos,xlims, ypos, ylims, zpos, zlims)


f2 = figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8]);
pcolor(x,t,u.'), shading interp, colormap(flipud(gray))
% title('Original Data','fontsize',18)
set_flat_figs(xlims, ylims, xpos_flat, ypos_flat)

f3 = figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8]);
% plot the edge/ridge-found points
pcolor(x,t,u.'), shading interp, colormap(flipud(gray))
hold on
scatter(xpts,tpts,'c','filled')
% title('Ridge Detection','fontsize',18)
set_flat_figs(xlims, ylims, xpos_flat, ypos_flat)


% plot initial spectral clustering
f4 = figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8]);
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
set_flat_figs(xlims, ylims, xpos_flat, ypos_flat)


f5 = figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8]);
% plot initial models
shifts = L*Asave{1};
pcolor(x,t,u.'), shading interp, colormap(flipud(gray))
hold on
for jj = 1:n
    if pct(jj) > 0.01
        xtmp = shifts(:,jj);
        plot(xtmp,t,'color',colors_mat{jj},'LineWidth',6)
    end
end
% title('Model Discovery','fontsize',18)
set_flat_figs(xlims, ylims, xpos_flat, ypos_flat)

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
set_flat_figs(xlims, ylims, xpos_flat, ypos_flat)


% plot final models
f7 = figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8]);
shifts = L*Asave{end};
pcolor(x,t,u.'), shading interp, colormap(flipud(gray))
hold on
for jj = 1:n
    if pct(jj) > 0.01
        xtmp = shifts(:,jj);
        plot(xtmp,t,'color',colors_mat{jj},'LineWidth',6)
    end
end
% title(['Final Models, Iteraton ',num2str(length(Asave))],'fontsize',18)
set_flat_figs(xlims, ylims, xpos_flat, ypos_flat)

%%

f8 = figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8]);
surfl(x,t,usrpca.'+upshift), shading interp, colormap(gray)
hold on
imagesc(x,t,usrpca.'/max(usrpca(:))), shading interp, colormap(flipud(gray))
set_3d_figs(views, xpos,xlims, ypos, ylims, zpos, zlims)
title('Shifted RPCA')

f9 = figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8]);
surfl(x,t,uspod.'+upshift), shading interp, colormap(gray)
hold on
imagesc(x,t,uspod.'/max(uspod(:))), shading interp, colormap(flipud(gray))
set_3d_figs(views, xpos,xlims, ypos, ylims, zpos, zlims)
title('Shifted POD')

f10 = figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8]);
surfl(x,t,upod.'+upshift), shading interp, colormap(gray)
hold on
imagesc(x,t,upod.'/max(upod(:))), shading interp, colormap(flipud(gray))
set_3d_figs(views, xpos,xlims, ypos, ylims, zpos, zlims)
title('Unshifted POD')

% print(f1,'figures/nls_data','-djpeg','-r600')
% print(f2,'figures/nls_data_flat','-djpeg','-r600')
% print(f3,'figures/nls_init_ridge','-djpeg','-r600')
% print(f4,'figures/nls_init_clusters','-djpeg','-r600')
% print(f5,'figures/nls_init_models','-djpeg','-r600')
% print(f6,'figures/nls_clusters','-djpeg','-r600')
% print(f7,'figures/nls_models','-djpeg','-r600')
% print(f8, 'figures/nls_sprca','-djpeg','-r600')
% print(f9, 'figures/nls_spod','-djpeg','-r600')
% print(f10, 'figures/nls_pod','-djpeg','-r600')

end



