function [RHS] = Prop_beta(t,Phi,A2inv,R,beta,N,dx)

nu = 10^(-3);

% % Calculate streamfunctions
 
u = Phi(1:N); u(1) = 0; 
v = Phi(N+1:end); v(1) = 0;

tau = A2inv*u; phi = A2inv*v;

Phi(1:N) = d_2x(N,dx)*tau; 
Phi(N+1:end) = d_2x(N,dx)*phi;

% Calculate w

RHS = 2*d_1x(N,dx)*Phi(N+1:end)+beta*d_1x(N,dx)*tau;

w = Omega_Solver(RHS,R,N,dx);

plot(w); pause(0.00001);

% Define Matrix

A = [zeros(N,N),-d_1x(N,dx); -d_1x(N,dx), zeros(N,N)];

b = -beta*[d_1x(N,dx)*tau; d_1x(N,dx)*phi];

c = [-w; zeros(N,1)];

d = -nu*[d_4x(N,dx),zeros(N,N);zeros(N,N),d_4x(N,dx)];

RHS = A*Phi+b+c; %+d*Phi;

end