%This program computes the bump functions c_k(1-x^2)^k, plots them and
%their power spectra. This is compared to a Gaussian cexp(-ax^2), which is
%scaled so that c*exp(-1)~10^(-14), Gaussian plots are red.
n= 8192
clear x y fy
x=linspace(-1,1,n+1);
a=12*log(10);%So sqrt(a/pi)*exp(-a)~10^(-14)
figure
hold on
for j=2:8
    y(1:n) = (1-x(1:n).*x(1:n)).^(2*j);
     s=2*sum(y)/n %Find the Riemann sum to normalize
    plot(x(1:n),y(1:n)/s)
end   
g(1:n)=exp(-a*x(1:n).*x(1:n));
  gs = 2*sum(g)/n
  g=g/gs;
  min(g)
plot(x(1:n),g(1:n),'Color',[1,0,0])
hold off
figure
hold on
for j=2:8
    y(1:n) = (1-x(1:n).*x(1:n)).^(2*j);
    fy = fft(y);
    s=2*sum(y)/n;
    plot(1:n,log10(abs(fy)/(n*s)))% plot Fourier power spectrum
end 
    fg = fft(g);
    plot(1:n,log10(abs(fg)/n),'Color',[1,0,0])% plot Fourier power spectrum
