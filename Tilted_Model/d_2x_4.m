function [A] = d_2x_4(N,dx)

% Calculate the second derivative of f - 4th accuracy

a = -5/(2*dx^2); a1 = 4/(3*dx^2); a2 = -1/(12*dx^2); 

dA = diag((a).*ones(1,N)); % diagonal entries
dAp1 = diag((a1).*ones(1,N-1),1); % i+1 entries
dAm1 = diag((a1).*ones(1,N-1),-1); % i-1 entries
dAp2 = diag((a2).*ones(1,N-2),2); % i+1 entries
dAm2 = diag((a2).*ones(1,N-2),-2); % i-1 entries

A = (dA + dAp1 + dAm1 + dAp2 + dAm2); 

A(1,end) = a1; A(1,end-1) = a2; A(2,end) = a2; %periodic BC
A(end,1) = a1; A(end,2) = a2; A(end-1,1) = a2; %periodic BC


end