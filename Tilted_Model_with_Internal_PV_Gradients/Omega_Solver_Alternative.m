function [w] = Omega_Solver_Alternative(RHS,R,N,dx)

% Define Matrix Operator

a = -2/dx^2-1; a1 = 1/dx^2; 

dA = diag((a).*ones(1,N)); % diagonal entries
dAp1 = diag((a1).*ones(1,N-1),1); % i+1 entries
dAm1 = diag((a1).*ones(1,N-1),-1); % i-1 entries

A = (dA + dAp1 + dAm1); 

A(1,end) = a1; %periodic BC
A(end,1) = a1; %periodic BC

% define tolerance

tol = 10^(-12);

% Initial Guess

w_old = zeros(N,1);
w_new = randn(N,1);
counter = 0;

while rms(w_new-w_old)>= tol
counter = counter+1;
% Define r

r = r_factor(w_new,R);

RHS_mod = RHS + d_2x(N,dx)*((1-r).*w_new);

% Invert

w_old = w_new;
w_new = A\RHS_mod;

end


w = w_new;
counter;
end