clear all; close all; clc

% Burgers' equation
%  u_t + u*u_x - eps*u_xx =0 

% setup
dt=0.1;
t=0:dt:20; eps=0.1; 
L=16; n=256; 
x2=linspace(-L/2,L/2,n+1); x=x2(1:n); k=(2*pi/L)*[0:n/2-1 -n/2:-1].';
k2=fftshift(k);

% initial data
u0=exp(-(x+2).^2).';
ut=fft(u0);  

[t,utsol]=ode45('burgers_rhs',t,ut,[],k,eps);
  
for j=1:length(t)
    u(:,j)=ifft( utsol(j,1:n).' );
end

figure(1), 
waterfall(x,t,real(u.')); colormap([0 0 0]);
view(15,35), set(gca,'Fontsize',[12])
set(gca,'Fontsize',[12],'Xlim',[-L/2 L/2],'Xtick',[-L/2 0 L/2],'Ylim',[0 10],'Ytick',[0 10],'Zlim',[0 1],'Ztick',[0 0.5 1])

u = abs(u);
x = x.';

folder = cd;
[folder,~,~] = fileparts(folder);
fname = fullfile(folder, 'burgers.mat');

save(fname, 'x', 't', 'u')


