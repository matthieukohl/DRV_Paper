function [A] = d_1x_f(N,dx)

% Calculate the second derivative of f

a = -3/(2*dx); a1p = 2/(dx); a2p = -1/(2*dx);


dA = diag((a).*ones(1,N)); % diagonal entries
dAp1 = diag((a1p).*ones(1,N-1),1); % i+1 entries
dAp2 = diag((a2p).*ones(1,N-2),2); % i-1 entries

A = (dA + dAp1 + dAp2); 

%periodic BC
A(end,1) = a1p; A(end,2) = a2p;
A(end-1,1) = a2p;

end