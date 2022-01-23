function F = root_solver_2(x,r)

% x(1) = sigma, x(2) = b, 

k1 = 1/sqrt(2*r)*sqrt(1-r*(x(1)^2+2)+sqrt(1+r^2*(x(1)^4+4*x(1)^2)-6*r*x(1)^2));

k2 = 1/sqrt(2*r)*sqrt(1-r*(x(1)^2+2)-sqrt(1+r^2*(x(1)^4+4*x(1)^2)-6*r*x(1)^2));

F(1) = tan(k2*x(2))-r*k1*k2/(1+x(1))*(k1/x(1)-1/(r*k1));

F(2) = tan(k1*x(2))-r*k1*k2/(1+x(1))*(k2/x(1)-1/(r*k2));


%F(1) = tan(x(2)/sqrt(r))+1/(sqrt(r)*(x(1)+1));
F(1) = x(2)-pi/2*sqrt(r)+(1+x(1))*r;
F(2) = x(2) - (1-x(1))/(x(1)*(1+x(1)));
%F(2) = x(2)-pi/2*sqrt(r)-(1+x(2))*r;
end