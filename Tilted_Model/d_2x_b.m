function [A] = d_2x_b(N,dx)

% Calculate the second derivative of f

a = 2/(dx^2); a1m = -5/(dx^2); a2m = 4/(dx^2); a3m = -1/dx^2;


dA = diag((a).*ones(1,N)); % diagonal entries
dAm1 = diag((a1m).*ones(1,N-1),-1); % i+1 entries
dAm2 = diag((a2m).*ones(1,N-2),-2); % i-1 entries
dAm3 = diag((a3m).*ones(1,N-3),-3); % i-1 entries

A = (dA + dAm1 + dAm2 + dAm3); 

%periodic BC
A(1,end) = a1m; A(1,end-1) = a2m; A(1,end-2) = a3m;
A(2,end) = a2m; A(2,end-1) = a3m;
A(3,end) = a3m;

end