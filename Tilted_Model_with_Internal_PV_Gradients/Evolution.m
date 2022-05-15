function [RHS] = Evolution(t,W,R,N,dx)

%plot(W(1:N)); pause(0.001);

% Calculate the r-factor

r = r_factor(W(1:N),R);

% Generate A and B Matrices

A = A_matrix(r,N,dx); B = B_matrix(r,N,dx); Binv = inv(B);

% Propagator matrix

M = [zeros(N,N),Binv;A,zeros(N,N)];

RHS = M*W;

end