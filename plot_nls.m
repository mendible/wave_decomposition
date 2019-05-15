function plot_for_paper3(params)
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
N = length(xpts);
colormap(flipud(gray))
colors_mat = {[0,0.4470,0.7410],...
              [0.8500,0.3250,0.0980],...
              [0.9290,0.6940,0.1250],'r','k','c'};

L = params.library(t);

f1=figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8])
% plot the 3D surface of the data
surfl(x,t,u.'+5), shading interp, colormap(gray)
hold on 
pcolor(x,t,u.'/max(u(:))), shading interp, colormap(flipud(gray))
view([12,23])
% title('Original Data','fontsize',18)
xlabel('x','fontsize',18)
ylabel('t','fontsize',18)
set(get(gca,'ylabel'),'rotation',0)
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'zticklabel',[])
set(f1,'Color','none')
print(f1,'figures/nls_data','-djpeg')


f2=figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8])
pcolor(x,t,u.'), shading interp, colormap(flipud(gray))
% title('Original Data','fontsize',18)
xlabel('x','fontsize',24)
ylabel('t','fontsize',24)
xlim([-10,10])
set(get(gca,'ylabel'),'rotation',0)
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'zticklabel',[])
set(f2,'Color',[0.5,0,0])
print(f2,'figures/nls_data_flat','-djpeg','-r600')

figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8])
% plot the edge/ridge-found points
pcolor(x,t,u.'), shading interp, colormap(flipud(gray))
hold on
scatter(xpts,tpts,'c','filled')
% title('Ridge Detection','fontsize',18)
xlabel('x','fontsize',24)
ylabel('t','fontsize',24)
xlim([-10,10])
set(get(gca,'ylabel'),'rotation',0)
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
box on
set(gcf,'Color','none')
print(gcf,'figures/nls_init_ridge','-djpeg')


% plot initial spectral clustering
figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8])
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
xlabel('x','fontsize',24)
ylabel('t','fontsize',24)
xlim([-10,10])
set(get(gca,'ylabel'),'rotation',0)
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
box on
set(gcf,'Color','none')
print(gcf,'figures/nls_init_clusters','-djpeg')


figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8])
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
xlabel('x','fontsize',24)
ylabel('t','fontsize',24)
xlim([-10,10])
set(get(gca,'ylabel'),'rotation',0)
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
box on
set(gcf,'Color','none')
print(gcf,'figures/nls_init_models','-djpeg')


% plot final spectral clustering
figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8])
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
xlabel('x','fontsize',24)
ylabel('t','fontsize',24)
xlim([-10,10])
set(get(gca,'ylabel'),'rotation',0)
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
box on
set(gcf,'Color',[0.5,0,0])
print(gcf,'figures/nls_clusters','-djpeg')


% plot final models
figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8])
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
xlabel('x','fontsize',24)
ylabel('t','fontsize',24)
xlim([-10,10])
set(get(gca,'ylabel'),'rotation',0)
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
box on
set(gcf,'Color','none')
print(gcf,'figures/nls_models','-djpeg')
end