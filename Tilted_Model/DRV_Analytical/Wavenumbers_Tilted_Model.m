% Tilted Model wavenumbers

sigma = 1;
x = linspace(0,1,100);

r = linspace(0.01,0.1,100);

k1 = 1./sqrt(2*r).*sqrt(-sqrt((-r*sigma^2-2*r+1).^2-4*r*sigma^2)-1+sigma^2*r+2*r);

k1approx = 1i./sqrt(r)-1i*(sigma^2+1)*sqrt(r)-1i/2*r.^(3/2)*(3*sigma^4+4*sigma^2+1);

k2 = 1./sqrt(2*r).*sqrt(sqrt((-r*sigma^2-2*r+1).^2-4*r*sigma^2)-1+sigma^2*r+2*r);

k2approx = 1i*sigma+1i*sigma*(sigma^2+1)*r-1i/2*sigma*(sigma^2+1)^2*r.^2;

figure(1)
plot(r,imag(k1)); hold on;
plot(r,imag(k1approx));
xlabel('r'); ylabel('Wavenumber k1')

figure(2)
plot(r,imag(k2)); hold on;
plot(r,imag(k2approx))
xlabel('r'); ylabel('Wavenumber k2')

