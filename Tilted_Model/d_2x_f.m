function [A] = d_2x_f(N,dx)

% Calculate the second derivative of f

a = 2/(dx^2); a1p = -5/(dx^2); a2p = 4/(dx^2); a3p = -1/dx^2;


dA = diag((a).*ones(1,N)); % diagonal entries
dAp1 = diag((a1p).*ones(1,N-1),1); % i+1 entries
dAp2 = diag((a2p).*ones(1,N-2),2); % i-1 entries
dAp3 = diag((a3p).*ones(1,N-3),3); % i-1 entries

A = (dA + dAp1 + dAp2 + dAp3); 

%periodic BC
A(end,1) = a1p; A(end,2) = a2p; A(end,3) = a3p;
A(end-1,1) = a2p; A(end-1,2) = a3p;
A(end-2,1) = a3p;

end