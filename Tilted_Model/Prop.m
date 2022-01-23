function [RHS] = Prop(t,Phi,A2inv,R,h,N,dx)
%L = 8*pi;
%x = linspace(0,L,N+1);
% % Calculate streamfunctions
 
u = Phi(1:N); u(1) = 0; 
v = Phi(N+1:end); v(1) = 0;

tau = A2inv*u; phi = A2inv*v;

Phi(1:N) = d_2x(N,dx)*tau; 
Phi(N+1:end) = d_2x(N,dx)*phi;

% Calculate w

RHS = 2*d_1x(N,dx)*Phi(N+1:end)-h*d_1x(N,dx)*phi;

w = Omega_Solver(RHS,R,N,dx);

% plot(w/norm(w)); xlabel('x'); ylabel('w');
% pause(0.0001)
% r = r_factor(w,R);
% PV2 = Phi(N+1:end) - Phi(1:N) + tau - mean(r.*w);

% plot(x(1:end-1),PV2); xlim([x(1) x(end-1)]); xlabel('x'); ylabel('PV');
% pause(0.0001)

% Define Matrix

A = [zeros(N,N),-d_1x(N,dx); -d_1x(N,dx), zeros(N,N)];

b = [d_1x(N,dx)*phi; d_1x(N,dx)*tau];

c = [-w; zeros(N,1)];

RHS = A*Phi+h*b+c; 

end