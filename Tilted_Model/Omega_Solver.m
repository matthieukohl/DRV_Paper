function [w] = Omega_Solver(RHS,R,N,dx)

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

% Define Matrix

M = Omega(r,N,dx);

% Invert

w_old = w_new;
w_new = M\RHS;

end


w = w_new;
counter;
end