function [RHS] = Prop_beta_1(t,Phi,R,beta,nu,N,dx)

% Calculate w

d = 2*d_3x(N,dx)*Phi(N+1:end)+beta*d_1x(N,dx)*Phi(1:N);

w = Omega_Solver(d,R,N,dx);

%plot(w); pause(0.01)

% Define Matrix

A = [-beta*d_1x(N,dx),-d_3x(N,dx); -d_3x(N,dx), -beta*d_1x(N,dx)];

c = [-w; zeros(N,1)];

d = -nu*[d_4x(N,dx),zeros(N,N);zeros(N,N),d_4x(N,dx)];

RHS = A*Phi+c+d*Phi;
RHS(1) = 0; RHS(N+1) = 0;

end