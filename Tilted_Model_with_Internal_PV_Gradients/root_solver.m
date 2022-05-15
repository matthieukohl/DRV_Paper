function F = root_solver(x,r,q)

% x(1) = sigma, x(2) = b, x(3) = c1, x(4) = c2

k1 = 1/sqrt(2*r)*sqrt(r*x(1)^2+sqrt((r*(x(1)^2+2)-1)^2-4*r*x(1)^2)+2*r-1);

k2 = 1/sqrt(2*r)*sqrt(r*x(1)^2-sqrt((r*(x(1)^2+2)-1)^2-4*r*x(1)^2)+2*r-1);

F(1) = x(3)*cosh(k1*x(2)/2)+x(4)*cosh(k2*x(2)/2);

F(2) = x(3)*k1*sinh(k1*x(2)/2)+x(4)*k2*sinh(k2*x(2)/2)+q*x(1)/(r*(x(1)+1));

F(3) = x(3)/k1 * sinh(k1*x(2)/2)+x(4)/k2*sinh(k2*x(2)/2)-q/(x(1)+1);

F(4) = 1/(k1^2-x(1)^2)*(-x(3)*k1*(r+1)*sinh(k1*x(2)/2)*(exp(-x(1)*x(2))+1)+...
    x(3)*x(1)*(r-1)*cosh(k1*x(2)/2)*(-exp(-x(1)*x(2))+1))+ ...
1/(k2^2-x(1)^2)*(-x(4)*k2*(r-1)*sinh(k2*x(2)/2)*(exp(-x(1)*x(2))+1) + ...
    x(4)*x(1)*(r-1)*cosh(k2*x(2)/2)*(exp(-x(1)*x(2))+1)) - q;

end
