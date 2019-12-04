function rhs=burgers_rhs(t,u0t,dummy,k,eps)

u0=ifft(u0t); u0x=ifft(i*k.*u0t);
rhs =-eps*(k.^2).*u0t - fft(u0x.*u0);
