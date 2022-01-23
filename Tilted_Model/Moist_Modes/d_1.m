function F = d_1(c0,L2)

% x(1) = c1, x(2) = L2

a20 = 1/sqrt(2)*sqrt(1-c0^2-sqrt(c0^4-6*c0^2+1));

b20 = 1/sqrt(2)*sqrt(1-c0^2+sqrt(c0^4-6*c0^2+1));

% Full System

F = (a20^2-b20^2)/(a20*tan(a20*L2)-b20*tan(b20*L2));



end