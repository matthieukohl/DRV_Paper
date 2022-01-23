clear; close all; 
sigma = linspace(0,1.62,100); 
r = linspace(0,1,50); 

[S,R] = meshgrid(sigma,r);
S = S'; R = R';

[k1,k2] = Wavenumbers(S,R);

RHS1 = -R.*k1.*k2./(S+1).*(1./(R.*k2)-S.*k2./(S.^2+R-1));

RHS2 = R.*k1.*k2./(1+S).*(-1./(R.*k1)+S.*k1./(S.^2+R-1));

figure(1)
contourf(R,S,imag(k1)); colorbar
xlabel('r'); ylabel('\sigma')
title('k1')

figure(2)
contourf(R,S,imag(k2)); colorbar
xlabel('r'); ylabel('\sigma')
title('k2')

figure(3)
contourf(R,S,imag(RHS2)); colorbar