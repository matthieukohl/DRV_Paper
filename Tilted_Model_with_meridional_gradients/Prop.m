function [RHS] = Prop(t,Phi,A2inv,R,N,dx,h1,h2,beta,alpha1,alpha2,alpha3)

% % Calculate streamfunctions
 
u = Phi(1:N); u(1) = 0; 
v = Phi(N+1:end); v(1) = 0;

tau = A2inv*u; phi = A2inv*v;

Phi(1:N) = d_2x(N,dx)*tau; 
Phi(N+1:end) = d_2x(N,dx)*phi;

% Calculate w

RHS = 2*d_1x(N,dx)*Phi(N+1:end)-1/2*(h1+h2)*d_1x(N,dx)*phi + 1/2*(h2-h1)*d_1x(N,dx)*tau...
+beta*d_1x(N,dx)*tau-(alpha3-alpha2)*d_2x(N,dx)*tau;

w = Omega_Solver(RHS,R,N,dx);

%w = Omega_Solver_galerkin(RHS,phi,R,N,dx);

%w = Omega_Solver_Alternative_Galerkin(RHS,phi,R,N,dx);

% plot(w/max(w));title(['h1=',num2str(h1),' ','h2=',num2str(h2),' ','r=',num2str(R)])
% pause(0.0001)
% q2 = d_2x(N,dx)*(phi-tau)+tau; %q2 = q2 -d_3x(N,dx)*(phi-tau);
% q1 = d_2x(N,dx)*(phi+tau)-tau; %q1 = q1 + d_3x(N,dx)*(phi+tau);
% plot(q2,'Color','r'); hold on;
% plot(q1+3*max(q2),'Color','b');
% hold off;
% %ylim([-0.1 1]); 
% title(['beta=',num2str(round(beta,2)),'  ','r=',num2str(round(R,2))])
% pause(0.0000001)


% Define Matrix

A = [zeros(N,N),-d_1x(N,dx); -d_1x(N,dx), zeros(N,N)];

b = [1/2*(h1 + h2)*d_1x(N,dx)*phi-1/2*(h2-h1)*d_1x(N,dx)*tau; -1/2*(h2-h1)*d_1x(N,dx)*phi+1/2*(h1+h2)*d_1x(N,dx)*tau];

c = [-w; zeros(N,1)];

d = beta*[d_1x(N,dx)*tau; d_1x(N,dx)*phi];

e = [alpha2*d_2x(N,dx)*tau; alpha1*d_2x(N,dx)*phi];

F = [d_2x(N,dx), zeros(N,N); zeros(N,N), d_2x(N,dx)];

RHS = A*Phi+b+c-d-e; 

% % Galilean Frame where u1 = 2U and u2 = 0
% % % Calculate streamfunctions
%  
% u = Phi(1:N); u(1) = 0; 
% v = Phi(N+1:end); v(1) = 0;
% 
% tau = A2inv*u; phi = A2inv*v;
% 
% Phi(1:N) = d_2x(N,dx)*tau; 
% Phi(N+1:end) = d_2x(N,dx)*phi;
% 
% % Calculate w
% 
% RHS = 2*d_1x(N,dx)*Phi(N+1:end)-h*d_1x(N,dx)*phi+beta*d_1x(N,dx)*tau...
%     -(alpha3-alpha2)*d_2x(N,dx)*tau;
% 
% w = Omega_Solver(RHS,R,N,dx);
% 
% %w = Omega_Solver_galerkin(RHS,phi,R,N,dx);
% 
% %w = Omega_Solver_Alternative_Galerkin(RHS,phi,R,N,dx);
% 
% plot(w/max(w));title(['beta=',num2str(beta),'','r=',num2str(R)])
% pause(0.001)
% % q2 = d_2x(N,dx)*(phi-tau)+tau; %q2 = q2 + t*d_1x(N,dx)*q2;
% % q1 = d_2x(N,dx)*(phi+tau)-tau; %q1 = q1 + t*d_1x(N,dx)*q1;
% % plot(q2,'Color','r'); hold on;
% % plot(q1+3*max(q2),'Color','b');
% % hold off;
% % %ylim([-0.1 1]); 
% % title(['beta=',num2str(round(beta,2)),'  ','r=',num2str(round(R,2))])
% % pause(0.0000001)
% 
% 
% % Define Matrix
% 
% A = [-d_1x(N,dx),-d_1x(N,dx); -d_1x(N,dx), -d_1x(N,dx)];
% 
% b = [d_1x(N,dx)*phi; d_1x(N,dx)*tau];
% 
% c = [-w; zeros(N,1)];
% 
% d = beta*[d_1x(N,dx)*tau; d_1x(N,dx)*phi];
% 
% e = [alpha2*d_2x(N,dx)*tau; alpha1*d_2x(N,dx)*phi];
% 
% RHS = A*Phi+h*b+c-d-e; 
% 



end