function F = EFT_solver(x,L2)

% x(1) = sigma, x(2) = L2

a2 = 1/sqrt(2)*sqrt(1-x(1)^2-sqrt(1-6*x(1)^2+x(1)^4));

b2 = 1/sqrt(2)*sqrt(1-x(1)^2+sqrt(1-6*x(1)^2+x(1)^4));

% sigma at O(1)

F(1) = (1-a2^2)/a2*tan(a2*L2)-(1-b2^2)/b2*tan(b2*L2);

end