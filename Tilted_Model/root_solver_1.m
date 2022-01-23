function F = root_solver_1(x,r)

% x(1) = sigma, x(2) = b, x(3) = c1, x(4) = c2

k1 = 1/sqrt(2*r)*sqrt(1-r*(x(1)^2+2)+sqrt(1+r^2*(x(1)^4+4*x(1)^2)-6*r*x(1)^2));

k2 = 1/sqrt(2*r)*sqrt(1-r*(x(1)^2+2)-sqrt(1+r^2*(x(1)^4+4*x(1)^2)-6*r*x(1)^2));

F(1) = x(3)*cos(k1*x(2))+x(4)*cos(k2*x(2));

F(2) = x(3)*k1*sin(k1*x(2))+x(4)*k2*sin(k2*x(2))-(x(1)-1)/r;

F(3) = (x(3)/k1)*sin(k1*x(2))+(x(4)/k2)*sin(k2*x(2))-(x(1)-1)/(x(1));

F(4) = x(3)*k1^2*cos(k1*x(2))+x(4)*k2^2*cos(k2*x(2))-(1-x(1)^2)/r;

end