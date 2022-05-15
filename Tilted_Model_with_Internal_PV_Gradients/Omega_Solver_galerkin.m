function [w] = Omega_Solver_galerkin(RHS,phi,R,N,dx)

% solve Omega-Equation for Galerkin approximation of 3-D Atmosphere

% define tolerance

tol = 10^(-5);

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

% Re-define the RHS

RHS_new = RHS-d_2x(N,dx)*(r.*d_1x(N,dx)*phi);

% Invert

w_old = w_new;
w_new = M\RHS_new;

end


w = w_new;
counter;
end