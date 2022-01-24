function [A] = Omega(r,N,dx)

% Initial Guess

% Define the Coefficients

a = -2*r/dx^2-1; a1 = r/dx^2; 

% Define the A-matrix 
    
dA = diag((a').*ones(1,N)); % diagonal entries
dAp1 = diag((a1(2:end)').*ones(1,N-1),1); % i+1 entries
dAm1 = diag((a1(1:end-1)').*ones(1,N-1),-1); % i-1 entries

A = (dA + dAp1 + dAm1);

% periodic BC

A(1,end) = a1(end);

A(end,1) = a1(1); 

end