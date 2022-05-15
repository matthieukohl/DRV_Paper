function [A] = d_4x_b(N,dx)

% Calculate the second derivative of f

a = 3/(dx^4); a1m = -14/(dx^4); a2m = 26/(dx^4); a3m = -24/dx^4;

a4m = 11/(dx^4); a5m = -2/(dx^4);


dA = diag((a).*ones(1,N)); % diagonal entries
dAm1 = diag((a1m).*ones(1,N-1),-1); % i+1 entries
dAm2 = diag((a2m).*ones(1,N-2),-2); % i-1 entries
dAm3 = diag((a3m).*ones(1,N-3),-3); % i-1 entries
dAm4 = diag((a4m).*ones(1,N-4),-4); % i-1 entries
dAm5 = diag((a5m).*ones(1,N-5),-5); % i-1 entries

A = (dA + dAm1 + dAm2 + dAm3 + dAm4 + dAm5); 

%periodic BC
A(1,end) = a1m; A(1,end-1) = a2m; A(1,end-2) = a3m; A(1,end-3) = a4m;A(1,end-4) = a5m;
A(2,end) = a2m; A(2,end-1) = a3m; A(2,end-2) = a4m; A(2,end-3) = a5m;
A(3,end) = a3m; A(3,end-1) = a4m; A(3,end-2) = a5m;
A(4,end) = a4m; A(4,end-1) = a5m;
A(5,end) = a5m;

end