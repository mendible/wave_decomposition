blue = [1, 123, 118]/255;
orange = [255, 82, 0]/255;
yellow = [0.9290,0.6940,0.1250];

err_pod = abs(u-params.ROM.upod);
err_spod = abs(u-params.ROM.uspod);
err_rpca = abs(u-params.ROM.usrpca);

top = max( [max(err_pod(:)), max(err_spod(:)), max(err_rpca(:)) ]);

figure

subplot(1,3,1)
pcolor(err_pod.'), shading flat, caxis([0 top]);
subplot(1,3,2)
pcolor(err_spod.'), shading flat, caxis([0 top]);
subplot(1,3,3)
pcolor(err_rpca.'), shading flat, caxis([0 top]);
colorbar;


figure
b = bar([1;2;3],[norm(err_pod); norm(err_spod); norm(err_rpca)]);
% set(gca,'xticklabel',{'POD','shifted POD','shifted RPCA'})
b.FaceColor = 'flat';

b.CData(1,:) = blue;
b.CData(2,:) = orange;
b.CData(3,:) = yellow;
set(gca,'fontsize',20)
yl = ylabel('$$l_2$$-norm error','Interpreter','latex');
pos = yl.Position;
pos(1) = -0.3;
set(yl,'position',pos)
xticks('')
