function [k1,k2] = Wavenumbers(s,r)

% Calculated the Wavenumbers in the Ascending Branch:
% w = c1*cos(k1*x)+c2*cos(k2*x)

inner_root = sqrt(1+r.^2.*(s.^4+4*s.^2)-6*r.*s.^2);

k1 = 1./sqrt(2*r).*sqrt(1-r.*(2+s.^2)+inner_root);
 
k2 = 1./sqrt(2*r).*sqrt(1-r.*(2+s.^2)-inner_root);
end

