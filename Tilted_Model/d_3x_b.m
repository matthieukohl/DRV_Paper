function [A] = d_3x_b(N,dx)

% Calculate the second derivative of f

a = 5/(2*dx^3); a1m = -9/(dx^3); a2m = 12/(dx^3); a3m = -7/dx^3;

a4m = 3/(2*dx^3);


dA = diag((a).*ones(1,N)); % diagonal entries
dAm1 = diag((a1m).*ones(1,N-1),-1); % i+1 entries
dAm2 = diag((a2m).*ones(1,N-2),-2); % i-1 entries
dAm3 = diag((a3m).*ones(1,N-3),-3); % i-1 entries
dAm4 = diag((a4m).*ones(1,N-4),-4); % i-1 entries

A = (dA + dAm1 + dAm2 + dAm3 + dAm4); 

%periodic BC
A(1,end) = a1m; A(1,end-1) = a2m; A(1,end-2) = a3m; A(1,end-3) = a4m;
A(2,end) = a2m; A(2,end-1) = a3m; A(2,end-2) = a4m;
A(3,end) = a3m; A(3,end-1) = a4m;
A(4,end) = a4m;

end