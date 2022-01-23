function F = EFT_solver_full(x,L2,r)

% x(1) = sigma, x(2) = L1

a1 = 1/sqrt(2*r)*sqrt(1-r*x(1)^2-sqrt(1-6*r*x(1)^2+r^2*x(1)^4));

b1 = 1/sqrt(2*r)*sqrt(1-r*x(1)^2+sqrt(1-6*r*x(1)^2+r^2*x(1)^4));

a2 = 1/sqrt(2)*sqrt(1-x(1)^2-sqrt(1-6*x(1)^2+x(1)^4));

b2 = 1/sqrt(2)*sqrt(1-x(1)^2+sqrt(1-6*x(1)^2+x(1)^4));


% a1 = 1/sqrt(2*r)*sqrt(r*x(1)^2-1+sqrt(1-6*r*x(1)^2+r^2*x(1)^4));
% 
% b1 = 1/sqrt(2*r)*sqrt(r*x(1)^2-1-sqrt(1-6*x(1)^2+r^2*x(1)^4));
% 
% a2 = 1/sqrt(2)*sqrt(x(1)^2-1+sqrt(1-6*x(1)^2+x(1)^4));
% 
% b2 = 1/sqrt(2)*sqrt(x(1)^2-1-sqrt(1-6*x(1)^2+x(1)^4));


% Full System

F(1) = r*b1*(b1^2-a1^2)*(tan(b2*L2)-b2/a2*tan(a2*L2))...
    +b2*(b2^2-a2^2)*(tan(b1*x(2))-b1/a1*tan(a1*x(2)));

F(2) = b1*(b1^2-a1^2)*(b2^2*tan(b2*L2)-a2*b2*tan(a2*L2))...
    +b2*(b2^2-a2^2)*(b1^2*tan(b1*x(2))-a1*b1*tan(a1*x(2)));

% F(1) = r*b1*(b1^2-a1^2)*(tanh(b2*L2)-b2/a2*tanh(a2*L2))...
%     +b2*(b2^2-a2^2)*(tanh(b1*x(2))-b1/a1*tanh(a1*x(2)));
% 
% F(2) = b1*(b1^2-a1^2)*(b2^2*tanh(b2*L2)-a2*b2*tanh(a2*L2))...
%     +b2*(b2^2-a2^2)*(b1^2*tanh(b1*x(2))-a1*b1*tanh(a1*x(2)));



end