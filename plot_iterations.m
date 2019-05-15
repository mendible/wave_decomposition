function plot_iterations(params,opt)

% opt.Ab = 'A' or 'B'
% opt.frame = frames to plot
% opt.sepFigs = if you want the same figure or separate figures
% opt.fname = video file name if desired

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
colors_mat = {'b','g','r','k','c'};
L = params.library(t);
count = 1;

colormap(flipud(gray))

for i = opt.frames
    if opt.sepFigs
        figure()
    else 
        figure(1)
        if count == 1
            figtitle = sgtitle('');
        end
    end
    figtitle.String = ['Iteration ',num2str(i)];
    pct = sum(Wsave{i})/sum(Wsave{i}(:));
    subplot(1,2,1)
    pcolor(x,t,u.'), shading interp, colormap(flipud(gray))
    hold on
    for jj = 1:n
        if pct(jj) > 0.01
            for kk = 1:N
                if Wsave{i}(kk,jj) == 1 
                    plot(xpts(kk),tpts(kk),'o','color',colors_mat{jj})
                end
            end
        end
    end
    
    if strcmp(opt.AB,'A')
        shifts = L*Asave{i};
    elseif strcmp(opt.AB,'B')
        shifts = L*Bsave{i};
    end
    
    subplot(1,2,2)
    pcolor(x,t,u.'), shading interp, colormap(flipud(gray))
    hold on
    for jj = 1:n
        if pct(jj) > 0.01
            xtmp = shifts(:,jj);
            plot(xtmp,t,colors_mat{jj},'LineWidth',4)
        end
    end
    hold off
    set(gcf,'position',[200 380 1040 400])
    drawnow
    if isfield(opt,'fname')
        F(count) = getframe(gcf) ;
        count = count+1;
    end
    figtitle.String = ('');
end

if isfield(opt,'fname')
    disp('saving video...')
    % create the video writer with 1 fps
    writerObj = VideoWriter(opt.fname, 'MPEG-4');
    writerObj.FrameRate = 10;
    % set the seconds per image
    % open the video writer
    open(writerObj);
    % write the frames to the video
    for j=1:length(F)
        % convert the image to a frame
        frame = F(j) ;
        writeVideo(writerObj, frame);
    end
    % close the writer object
    close(writerObj);
end