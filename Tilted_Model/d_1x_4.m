function [A] = d_1x_4(N,dx)

% Calculate the second derivative of f

a = 0; a1p = 2/(3*dx); a1m = -a1p; a2p = -1/(12*dx); a2m = -a2p;

dA = diag((a).*ones(1,N)); % diagonal entries
dAp1 = diag((a1p).*ones(1,N-1),1); % i+1 entries
dAm1 = diag((a1m).*ones(1,N-1),-1); % i-1 entries
dAp2 = diag((a2p).*ones(1,N-2),2); % i+1 entries
dAm2 = diag((a2m).*ones(1,N-2),-2); % i-1 entries

A = (dA + dAp1 + dAm1 + dAp2 + dAm2); 

A(1,end) = a1m; A(1,end-1) = a2m; %periodic BC
A(2,end) = a2m;
A(end,1) = a1p; A(end,2) = a2p;%periodic BC
A(end-1,1) = a2p;

end