clear all; close all; clc

L=30; n=1024;
x2=linspace(-L/2,L/2,n+1); x=x2(1:n); 
k=(2*pi/L)*[0:n/2-1 -n/2:-1].';
t=linspace(0,2*pi,501);

omega=2;
u1=2*sech(x+7).*exp(i*omega*x);
u2=2*sech(x-7).*exp(-i*omega*x);
u=u1+u2;
ut=fft(u);

[t,utsol]=ode45('nls_rhs',t,ut,[],k);

for j=1:length(t)
   usol(j,:)=ifft(utsol(j,:));
end

figure(1), subplot(1,2,1),
surfl(x,t,abs(usol)); shading interp, colormap(gray);
title('original data')

u = abs(usol).';
x = x.';

save nls.mat x t u 