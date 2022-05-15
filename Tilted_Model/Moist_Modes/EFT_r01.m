% find the max growth rate and half-ascent
% for EFT at r=0.1
close all; clear;

r = 0.1;

N = 400;

L2 = linspace(1.3,3,N);
L2 = flip(L2);


x0 = [0.9,0.6];
for ii =1:N
if ii>2
x0(1) = 2*s(ii-1)-s(ii-2); 
x0(2) = 2*b(ii-1)-b(ii-2); 
end    
    
[sol,fval(ii,:)] = fsolve(@(x) EFT_solver_paper(x,L2(ii),r),x0);
s(ii) = sol(1); b(ii) = sol(2);
end

figure(1)
plot(L2,s); 
xlabel('L2');
ylabel('\sigma')

figure(2)
plot(L2,b);
xlabel('L2');
ylabel('b')

figure(3)
plot(L2,real(fval(:,1))); hold on;
plot(L2,real(fval(:,2)));

figure(4)
plot(L2,imag(fval(:,1))); hold on;
plot(L2,imag(fval(:,2)));

[pholder,ind] = max(s);
disp('max \sigma=')
s(ind)
disp('max b=')
b(ind)