function [k1,k2] = Wavenumbers_Complex(sigma,r)

% Calculate the Wavenumbers in the Ascending Branch:
% w = c1*cos(k1*x)+c2*cosh(k2*x)


%k1 = 1/sqrt(2*r)*sqrt(1-r*(2+sigma^2)+sqrt(1+r^2*(sigma^4+4*sigma^2)-6*r*sigma^2));


%k2 = 1/sqrt(2*r)*sqrt(r*(2+sigma^2)-1+sqrt(1+r^2*(sigma^4+4*sigma^2)-6*r*sigma^2));


k1 = 1./sqrt(2*r).*sqrt(1-r.*(2+sigma.^2)+sqrt(1+r.^2.*(sigma.^4+4*sigma.^2)-6*r.*sigma.^2));


k2 = 1./sqrt(2*r).*sqrt(r.*(2+sigma.^2)-1+sqrt(1+r.^2.*(sigma.^4+4*sigma.^2)-6*r.*sigma.^2));


end