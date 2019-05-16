% kdv.m - Solve KdV equation by Fourier spectral/ETDRK4 scheme
%         A.-K. Kassam and L. N. Trefethen 4/03
%
% https://www.mathworks.com/matlabcentral/answers/181591-code-to-solve-kdv-ecuation-with-an-animation-of-2-solitions
%
% This code solves the Korteweg-de Vries eq. u_t+uu_x+u_xxx=0
% with periodic BCs on [-pi,pi] and initial condition given by
% a pair of solitons.  The curve evolves up to t=0.005 and at
% the end u(x=0) is printed to 6-digit accuracy.  Changing N
% to 384 and h to 2.5e-7 improves this to 10 digits but takes
% four times longer.
% Set up grid and two-soliton initial data:
N = 512;
x = (2*pi/N)*(-N/2:N/2-1)';
tf = 0.6;
u0 = zeros(size(x));% may have to switch these
u = u0;
alpha = 25;
beta = 1;
kappa = 1;
speed = (alpha+2*beta)/3
coeff = sqrt((alpha-beta)/(12*kappa))
u0 = beta+(alpha-beta)*sech(coeff*x).^2

%     u = u + 3*A^2*sech(.5*(A*(x+rand*2))).^2;
% u = u + 3*A^2*sech(.5*(A*(x))).^2; %2/24
% % % beta = 1.5;
% % % alpha = 12*beta^2;
% % % speed = 4*beta^1;
% % % u = alpha*sech(beta*(x)).^2;

p = plot(x,u0,'linewidth',3);
axis([-pi pi -10 50]), grid on
% Precompute ETDRK4 scalar quantities (Kassam-Trefethen):
h = 1e-5;                               % time step
k = [0:N/2-1 0 -N/2+1:-1]';             % wave numbers
L = 1i*k.^3;                            % Fourier multipliers
E = exp(h*L); E2 = exp(h*L/2);
M = 64;                                 % no. pts for complex means
r = exp(2i*pi*((1:M)-0.5)/M);           % roots of unity
LR = h*L(:,ones(M,1))+r(ones(N,1),:);
Q  = h*mean(                  (exp(LR/2)-1)./LR   ,2);
f1 = h*mean((-4-LR+exp(LR).*(4-3*LR+LR.^2))./LR.^3,2);
f2 = h*mean(    (4+2*LR+exp(LR).*(-4+2*LR))./LR.^3,2);
f3 = h*mean((-4-3*LR-LR.^2+exp(LR).*(4-LR))./LR.^3,2);
g = -.5i*k;
% Time-stepping by ETDRK4 formula (Cox-Matthews):
set(gcf,'doublebuffer','on')
%disp('press <return> to begin'), pause  % wait for user input
t_stp = 0; step = 0; v = fft(u0); t = t_stp;
while t_stp+h/2 < tf
    step = step+1;
    t_stp = t_stp+h;
    %     tsave = [tsave,t];
    Nv = g.*fft(real(ifft(v)).^2);
    a = E2.*v+Q.*Nv;        Na = g.*fft(real(ifft(a)).^2);
    b = E2.*v+Q.*Na;        Nb = g.*fft(real(ifft(b)).^2);
    c = E2.*a+Q.*(2*Nb-Nv); Nc = g.*fft(real(ifft(c)).^2);
    v = E.*v+(Nv.*f1+(Na+Nb).*f2+Nc.*f3);
    u0 = real(ifft(v));
%     usave = [usave,u];
    if mod(step,100)==0
        t = [t,t_stp];
        u = [u,u0];
        set(p,'ydata',u0)
        title(sprintf('t = %7.5f',t_stp),'fontsize',18), drawnow
    end
end
% text(-2.4,900,sprintf('u(0) = %11.7f',u(N/2+1)),...
%     'fontsize',18,'color','r')

t = t.';

folder = cd;
[folder,~,~] = fileparts(folder);
fname = fullfile(folder, 'kdv.mat');

save(fname, 'x', 't', 'u')