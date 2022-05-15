function [Phi_final] = Implicit_Stepper(Phi0,R,N,dx,tN,dt)

nu = 0;

% Pre-define

Phi_final = zeros(2*N,tN);

% Define Operator

M1 = [zeros(N,N),-d_1x(N,dx);-d_1x(N,dx),zeros(N,N)];

M2 = [-nu*d_4x(N,dx),zeros(N,N);zeros(N,N),-nu*d_4x(N,dx)];

Operator = inv(diag(ones(1,2*N))-dt*(M1+M2));

% Initialize

Phi = Phi0;

for tt = 1:tN
tt
Phi = rescale(Phi);
% Calculate-w

RHS = 2*d_1x(N,dx)*Phi(N+1:end);
w = Omega_Solver(RHS,R,N,dx);
c = [-w;zeros(N,1)];

% Propagate Implicitly
    
Phi = Operator*(Phi+dt*c);
Phi_final(:,tt) = Phi;
    
end

end