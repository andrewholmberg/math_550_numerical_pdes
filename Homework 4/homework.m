a = -0.5;
b = 0.5;
N = 64;
L = b-a;
h = L/N;
x = a+h*(1:N)';

sig = 0.1;
f = @(x) (1/sig/sqrt(2*pi))*exp(-0.5*(x./sig).^2);
d2_f = @(x) (1/(sig^5*sqrt(2*pi)))*(x.^2-sig^2).*exp(-0.5*(x./sig).^2);


%%%%%%%%%%%
kk = [0:N/2-1 0 -N/2+1:-1]';
kk2 = [0:N/2 -N/2+1:-1]';
ik = ((2*pi)/L)*1i*kk;
ik2 = ((2*pi)/L)*1i*kk2;
%%%%%%%%%%%

d_xx_fft_1 = ifft(ik.^2.*fft(f(x)));
d_xx_fft_2 = ifft(ik2.^2.*fft(f(x)));
d_xx_exact = d2_f(x);           
plot(x,d_xx_fft_1,'-o','linewidth',3)
hold all
plot(x,d_xx_fft_2,'-p','linewidth',3)

hold all

plot(x,d_xx_exact,'--k','linewidth',3)