function [B] = d_3x_4(N,dx)

% Finite difference definition of d_3x operator with second accuracy

b = 0; b1p = -13/(8*dx^3); b1m = -b1p; b2p = 1/(dx^3); b2m = -b2p;

b3p = -1/(8*dx^3); b3m = -b3p;

dB = diag((b).*ones(1,N)); % diagonal entries
dBp1 = diag((b1p).*ones(1,N-1),1); % i+1 entries
dBm1 = diag((b1m).*ones(1,N-1),-1); % i-1 entries
dBp2 = diag((b2p).*ones(1,N-2),2); % i+2 entries
dBm2 = diag((b2m).*ones(1,N-2),-2); % i-2 entries
dBp3 = diag((b3p).*ones(1,N-3),3); % i+2 entries
dBm3 = diag((b3m).*ones(1,N-3),-3); % i-2 entries

B = (dB + dBp1 + dBm1+ dBp2 + dBm2 + dBp3 + dBm3);

B(1,end) = b1m; B(1,end-1) = b2m; B(1,end-2) = b3m;
B(2,end) = b2m; B(2,end-1) = b3m;
B(3,end) = b3m;

B(end,1) = b1p; B(end,2) = b2p; B(end,3) = b3p;
B(end-1,1) = b2p; B(end-1,2) = b3p;
B(end-2,1) = b3p;

end