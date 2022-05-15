function F = root_solver_beyond_cut_off(x,r)

% x(1) = sigma, x(2) = b

k1 = 1/sqrt(2*r)*sqrt(-(1-r*(x(1)^2+2)+sqrt(1+r^2*(x(1)^4+4*x(1)^2)-6*r*x(1)^2)));

k2 = 1/sqrt(2*r)*sqrt(-(1-r*(x(1)^2+2)-sqrt(1+r^2*(x(1)^4+4*x(1)^2)-6*r*x(1)^2)));



% % % % Solve the two tangent equations
% % % % 
F(1) = tanh(k1*x(2))-r*k1*k2/(1+x(1))*(-1/(r*k2)-k2*x(1)/(x(1)^2+r-1));
% 
F(2) = tanh(k2*x(2))-r*k1*k2/(1+x(1))*(-1/(r*k1)-x(1)*k1/(x(1)^2+r-1));


end