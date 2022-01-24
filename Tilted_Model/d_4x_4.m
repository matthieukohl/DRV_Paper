function [A] = d_4x_4(N,dx)

% Calculate the fourth-order derivative of f with 4th accuracy

a = 28/(3*dx^4); a1 = -13/(2*dx^4); a2 = 2/dx^4; a3 = -1/(6*dx^4);

dA = diag((a).*ones(1,N)); % diagonal entries
dAp1 = diag((a1).*ones(1,N-1),1); % i+1 entries
dAm1 = diag((a1).*ones(1,N-1),-1); % i-1 entries
dAp2 = diag((a2).*ones(1,N-2),2); % i+2 entries
dAm2 = diag((a2).*ones(1,N-2),-2); % i-2 entries
dAp3 = diag((a3).*ones(1,N-3),3); % i+2 entries
dAm3 = diag((a3).*ones(1,N-3),-3); % i-2 entries


A = (dA + dAp1 + dAm1 + dAp2 + dAm2 + dAp3 + dAm3); 

A(1,end) = a1; A(1,end-1) = a2; A(1,end-2) = a3; A(2,end) = a2; A(2,end-1)= a3; A(3,end) = a3; %periodic BC
A(end,1) = a1; A(end,2) = a2; A(end,3) = a3; A(end-1,1) = a2; A(end-1,2) = a3; A(end-2,1) = a3; % periodic BC

end