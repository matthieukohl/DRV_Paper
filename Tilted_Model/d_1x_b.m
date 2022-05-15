function [A] = d_1x_b(N,dx)

% Calculate the second derivative of f

a = 3/(2*dx); a1m = -2/(dx); a2m = 1/(2*dx);


dA = diag((a).*ones(1,N)); % diagonal entries
dAm1 = diag((a1m).*ones(1,N-1),-1); % i+1 entries
dAm2 = diag((a2m).*ones(1,N-2),-2); % i-1 entries

A = (dA + dAm1 + dAm2); 

%periodic BC
A(1,end) = a1m; A(1,end-1) = a2m;
A(2,end) = a2m;

end