function F = EFT_solver_1(x,c0,L2)

% x(1) = sigma, 

a20 = 1/sqrt(2)*sqrt(c0^2-2+sqrt(c0^4-12*c0^2+4));

a21 = 1/(4*a20)*(2*c0+2*c0*(c0^2-6))/sqrt(c0^4-12*c0^2+4)*x(1);

b20 = 1/sqrt(2)*sqrt(c0^2-2-sqrt(c0^4-12*c0^2+4));

b21 = -1/(4*b20)*(2*c0+2*c0*(c0^2-6))/sqrt(c0^4-12*c0^2+4)*x(1);

% Full System

F(1) = -2*(b20/a20*a21*L2*sech(L2*a20)^2+(b21/a20-a21*b20/a20^2)*tanh(a20*L2) ...
    -b21*L2*sech(b20*L2)^2)-(b20*a21+b21*a20)*tanh(a20*L2)-...
    b20*a20*a21*L2*sech(a20*L2)^2+2*b20*b21*tanh(b20*L2)...
    +b20^2*b21*L2*sech(b20*L2)^2+pi/(2*sqrt(2))*b20*(b20^2-a20^2);

end