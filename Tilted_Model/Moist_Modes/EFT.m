% Check EFT result for asymptotic sigma at O(1); equation C4
close all;
clear;
N = 100;
L2 = linspace(2,5,N);
sigma = zeros(N,1); 

x0 = [1.3];

for ii=1:N
sol = fsolve(@(x) EFT_solver(x,L2(ii)),x0);
sigma(ii) = sol(1); 
end

figure(1)
plot(L2,sigma); hold on;
xlabel('L2'); ylabel('sigma')

% Calculate O(sqrt(r)) results

s0 = max(sigma); L_2 = L2(find(sigma==max(sigma)));
x0 = 3;
sol = fsolve(@(x) EFT_max_solver(x,s0,L_2),x0);
s1 = sol(1);

R = linspace(0,0.01,N); s_theory = s0+sqrt(R)*s1;

figure(10)
semilogx(R,s_theory)
xlabel('r'); 
ylabel('sigma')


% Check the asypmtotic results

rr = linspace(0.00001,0.01,N);
sigma_num = zeros(1,N); L1_num = zeros(N,1);
sigma_theory = s0+sqrt(rr)*real(s1);
d1 = d_1(s0,L_2); d2 = d_2(s0,real(s1),d1,L_2);
L1_theory = pi/(2)*sqrt(rr)+d1*rr+d2*rr.^(3/2);

for ii =1:N
x0 = [sigma_theory(ii),L1_theory(ii)];
sol = fsolve(@(x) EFT_solver_full(x,L_2,rr(ii)),x0);
sigma_num(ii) = sol(1); L1_num(ii) = sol(2);
end

% Plot Results

figure(20)
plot(rr,L1_num); hold on;
plot(rr,L1_theory)
xlabel('r');
ylabel('b');
legend('Numerical','Theory','Location','Northwest')
legend boxoff
set(gca,'TickDir','out','Box','off','Layer','top')
set(gca,'fontsize',14)
set(gca,'linewidth',1.5)
saveas(gcf,'EFT_L1.png')

figure(21)
plot(rr,sigma_num); hold on;
plot(rr,sigma_theory)
xlabel('r');
ylabel('sigma');
legend('Numerical','Theory')
legend boxoff
set(gca,'TickDir','out','Box','off','Layer','top')
set(gca,'fontsize',14)
set(gca,'linewidth',1.5)
saveas(gcf,'EFT_sigma.png')


% Comparison DRV vs. Moist Modes

R = linspace(0,1,100);
s_DRV = 1.6180-2.9756*sqrt(R);
b_DRV = pi/2*sqrt(R)+2.6180*R;

s_MM = 1.0494-2.1272*sqrt(R);
b_MM = pi/2*sqrt(R)+1.4837*R;

figure(30)
plot(R,b_DRV); hold on;
plot(R,b_MM)
xlabel('r');
ylabel('b');
legend('DRV','Moist Modes')
legend boxoff
set(gca,'TickDir','out','Box','off','Layer','top')
set(gca,'fontsize',14)
set(gca,'linewidth',1.5)
saveas(gcf,'DRV_MM_L1.png')

figure(31)
plot(R,s_DRV); hold on;
plot(R,s_MM)
xlabel('r');
ylabel('sigma');
legend('DRV','Moist Modes')
legend boxoff
set(gca,'TickDir','out','Box','off','Layer','top')
set(gca,'fontsize',14)
set(gca,'linewidth',1.5)
saveas(gcf,'DRV_MM_sigma.png')

