function [A] = d_4x_f(N,dx)

% Calculate the second derivative of f

a = 3/(dx^4); a1p = -14/(dx^4); a2p = 26/(dx^4); a3p = -24/dx^4;

a4p = 11/(dx^4); a5p = -2/(dx^4);


dA = diag((a).*ones(1,N)); % diagonal entries
dAp1 = diag((a1p).*ones(1,N-1),1); % i+1 entries
dAp2 = diag((a2p).*ones(1,N-2),2); % i-1 entries
dAp3 = diag((a3p).*ones(1,N-3),3); % i-1 entries
dAp4 = diag((a4p).*ones(1,N-4),4); % i-1 entries
dAp5 = diag((a5p).*ones(1,N-5),5); % i-1 entries

A = (dA + dAp1 + dAp2 + dAp3 + dAp4 + dAp5); 

%periodic BC
A(end,1) = a1p; A(end,2) = a2p; A(end,3) = a3p; A(end,4) = a4p; A(end,5) = a5p;
A(end-1,1) = a2p; A(end-1,2) = a3p; A(end-1,3) = a4p; A(end-1,4) = a5p;
A(end-2,1) = a3p; A(end-2,2) = a4p; A(end-2,3) = a5p;
A(end-3,1) = a4p; A(end-3,2) = a5p;
A(end-4,1) = a5p;

end