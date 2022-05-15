function [A] = d_1x(N,dx)

% Calculate the second derivative of f

a = 0; a1p = 1/(2*dx); a1m = -1/(2*dx);

dA = diag((a).*ones(1,N)); % diagonal entries
dAp1 = diag((a1p).*ones(1,N-1),1); % i+1 entries
dAm1 = diag((a1m).*ones(1,N-1),-1); % i-1 entries

A = (dA + dAp1 + dAm1); 

A(1,end) = a1m; %periodic BC
A(end,1) = a1p; %periodic BC

end