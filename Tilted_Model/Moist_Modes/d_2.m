function F = d_2(c0,c1,d1,L2)

% x(1) = c1, x(2) = L2

a20 = 1/sqrt(2)*sqrt(1-c0^2-sqrt(c0^4-6*c0^2+1));

a21 = 1/(4*a20)*(-2*c0-2*c0*(c0^2-3)/sqrt(c0^4-6*c0^2+1))*c1;

b20 = 1/sqrt(2)*sqrt(1-c0^2+sqrt(c0^4-6*c0^2+1));

b21 = 1/(4*b20)*(-2*c0+2*c0*(c0^2-3)/sqrt(c0^4-6*c0^2+1))*c1;

% Full System

F = pi/2*c0^2+d1^2*(a21*tan(a20*L2)+a20*a21*L2*sec(a20*L2)^2 ...
    -b21*tan(b20*L2)-b20*b21*L2*sec(b20*L2)^2)/(b20^2-a20^2);



end