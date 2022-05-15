function [A] = d_1xs(N,x)

% Calculate the first derivative of f

dxi = circshift(diff(x),1); dxip1 = diff(x); 

a = 0; a1p = 1./(dxi+dxip1); a1m = -a1p;

dA = diag((a).*ones(1,N)); % diagonal entries
dAp1 = diag((a1p(1:N-1)).*ones(1,N-1),1); % i+1 entries
dAm1 = diag((a1m(2:N)).*ones(1,N-1),-1); % i-1 entries

A = (dA + dAp1 + dAm1); 

A(1,end) = a1m(1); %periodic BC
A(end,1) = a1p(end); %periodic BC

end