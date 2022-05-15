% Solve the Mode-Equation as an initial value problem:
% d_tt ((rw)_xx-4*w) = (rw)_xxxx+4*w_xx
% d_tt B*w = A*w -> d_tt w = B^(-1)*A*w

clear;
close all;

% Define Grid

N = 400; 
L = 16*pi; 
x = linspace(0,L,N+1);
dx = L/(N);

% Define time-step

t_final = 400;
tN = 16000;
dt = t_final/tN;

% define r-factor

R = [1.0,0.8,0.6,0.4,0.3,0.2,0.1,0.05,0.02,0.01];

% define beta-factor

beta = 0;

% define hyperviscosity

nu = 0; %1*10^(-3);

% define Initial Conditions

IC = 'rand';
%IC = 'sinusoidal';

% define variables

sigma = zeros(length(R),1); lambda = zeros(length(R),1); 

w_final = zeros(N+1,length(R));

%[partial_t,adv_up,adv_low,adv_bp_up,adv_bp_low,adv_pb_up,adv_pb_low,heating,...
%budget_w,budget_w_mode,budget_PV] = deal(zeros(length(R),1));

[partial_t_amp,bp_amp,pb_amp,heat_amp] = deal(zeros(length(R),1));

[w_equation,w_equation_mode,PV_equation,...
phi_equation,tau_equation] = deal(zeros(N,length(R))); 

% % pre-define differentiation operators
% 
% % 2nd-order accuracy
% 
% d_1 = d_1x(N,dx); d_2 = d_2x(N,dx); d_3 = d_3x(N,dx); d_4 = d_4x(N,dx);
% 
% % 4th-order accuracy
% 
% d_1 = d_1x_4(N,dx); d_2 = d_2x_4(N,dx); d_3 = d_3x_4(N,dx); d_4 = d_4x_4(N,dx);

% define the mass matrix

M = d_2x(N,dx); M(1,:) = 1;

LHS = [M zeros(N,N); zeros(N,N) M];

% specify ODE-Values

options = odeset('Mass', LHS,'Events',@myEvent);

%options = odeset('Mass', LHS);

for ii = 1:1
ii

% Initialize Values

[w_test,r_test,tau_test,phi_test] = deal(zeros(N,3));

[Phi_test] = deal(zeros(2*N,3));

% Initial Values for Phi

if strcmp(IC,'rand') == 1
% random-guess
phi =   randn(1,N); phi = phi - mean(phi); 
tau =   randn(1,N); tau = tau - mean(tau); 
else

% sinusoidal-guess

% phi = sin(x)+sin(3*x)+sin(5*x); phi = phi(1:end-1);
% tau = sin(x)+sin(3*x)+sin(5*x); tau = tau(1:end-1);

phi = sin(x)+sin(2*x)+sin(4*x)+sin(6*x)+sin(8*x)+sin(10*x); phi = phi(1:end-1);
tau = sin(x)+sin(2*x)+sin(4*x)+sin(6*x)+sin(8*x)+sin(10*x); tau = tau(1:end-1);

end

Phi0 = [tau,phi]';

% Define inverse matrix

A2 = d_2x(N,dx); A2(1,:) = 1; A2inv = inv(A2);

%Propagator 
total_time = 0;
g = [];
tt = [];
tte = [];

%[t,y] = ode45(@(t,y) Prop_beta_1(t,y,A2inv,R(ii),beta,N,dx),[0 t_final],Phi0,options);

while total_time <= t_final
[t,y,te,ye,ie] = ode45(@(t,y) Prop_beta_1(t,y,R(ii),beta,nu,N,dx),[0 t_final],Phi0,options);
Phi0 = y(end,:)/100;
total_time = t(end)+total_time;
g = cat(1,g,y);
tt = cat(1,tt,t);
tte = cat(1,tte,te);
end

dt = t(end)-t(end-1);

% Retrieve phi/tau fields

phi = y(:,N+1:end)'; tau = y(:,1:N)';

phi_end = phi(:,end); phi_p_end = phi(:,end-1);
tau_end = tau(:,end); tau_p_end = tau(:,end-1); 

% Last step w-inersion

RHS = 2*d_3x(N,dx)*phi(:,end)+beta*d_1x(N,dx)*tau(:,end);
w = Omega_Solver(RHS,R(ii),N,dx);
r = r_factor(w,R(ii));

% final fields 

w_final(:,ii) = [w;w(1)];

% Calculate growthrate

sigma(ii) = growthrate(y,t);

% Calculate Asymmetry

lambda(ii) = asymmetry(w_final(:,ii));

% % check-equation closure
 
w_equation(:,ii) = 1/norm(w)*(d_2x(N,dx)*(r.*w)-w-RHS);

tau_equation(:,ii) = 1/norm(tau(:,end))*(1/dt*d_2x(N,dx)*(tau(:,end)-tau(:,end-1))+d_3x(N,dx)*phi(:,end)+beta*d_1x(N,dx)*tau_end+w);

phi_equation(:,ii) = 1/norm(phi(:,end))*(1/dt*d_2x(N,dx)*(phi(:,end)-phi(:,end-1))+d_3x(N,dx)*tau(:,end)+beta*d_1x(N,dx)*phi_end);

% Calculate the PV-Budget

%[partial_t_term,low,bp_low,pb_low,heat] = main_PV_beta_1_test(w,tau_end,tau_p_end,phi_end,phi_p_end,R(ii),sigma(ii),beta,dt,N,dx);

[partial_t_term,low,bp_low,pb_low,heat,diffusion] = main_PV_beta_1(w,tau_end,tau_p_end,phi_end,phi_p_end,R(ii),beta,nu,dt,N,dx);

PV_equation(:,ii) = (partial_t_term + pb_low + bp_low - heat);

% Calculate norm of PV-Budget terms

partial_t_amp(ii) = rms(partial_t_term)/max(bp_low);
bp_amp(ii) = rms(bp_low)/max(bp_low);
pb_amp(ii) = rms(pb_low)/max(bp_low);
heat_amp(ii) = rms(heat)/max(bp_low);


figure(24)
plot(x,-[heat;heat(1)]); hold on;
plot(x,[partial_t_term;partial_t_term(1)]); hold on
plot(x,[bp_low;bp_low(1)]); hold on;
plot(x,[pb_low;pb_low(1)]); hold on
plot(x,[PV_equation(:,ii);PV_equation(1,ii)])
plot(x,w_final(:,ii)); hold on;
xlabel('x'); ylabel('Magnitude')
legend('heat','Partial t','uqx','vqy','Residual','w');
%title(['r=',num2str(R(ii)),' ','beta=',num2str(beta),' ','dx=',num2str(round(dx,2)),' ','\sigma=',num2str(round(sigma(ii),2)),' ','nu=',num2str(nu),' ',IC])
%savefig(['/net/aimsir/archive1/mkohl/Mode_Zurita_Gotor/test/nu_53/PV_Budget_',IC,'_nu_r',num2str(R(ii)),'_b',num2str(beta),'_dx',num2str(round(dx,2)),'.fig'])
title(['r=',num2str(R(ii))])
%saveas(gcf,['Plots_Proposal/Large_Domain/PV_Budget_',IC,'_r',num2str(R(ii)),'_b',num2str(beta),'_dx',num2str(round(dx,2)),'.png'])
close all;
end

% Vertical Profile - r = 0.1

figure(1)
plot(x,w_final(:,7)/norm(w_final(:,7)),'linewidth',1.6); hold on
xlabel('x','interpreter','latex');
ylabel('w','interpreter','latex')
ax = gca;
ax.XTick = [0,2*pi,4*pi,6*pi,8*pi];
ax.XTickLabel = {'0','2\pi','4\pi','6\pi','8\pi'};
title('w-profile for r=0.1')
set(gca,'fontsize', 14);
xlim([0 8*pi])
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'linewidth',1.5)
%saveas(gcf,['Plots_Proposal/Large_Domain/w_0.1_',IC,'_dx',num2str(round(dx,2)),'.png'])


% Modal Theory vs. GCM-asymmetry

r_factor = [0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.6,0.8,1.0];

lambda_omega_mode = [0.9283,0.9173,0.8996,0.8313,0.7820,0.7343,0.6871,...
0.5674,0.5187,0.4994];

figure(2)
semilogx(r_factor,lambda_omega_mode,'Color','red','linewidth',1.6); hold on;
semilogx(R,lambda,'Color','blue','linewidth',1.6); hold on;
xlabel('Reduction factor r','interpreter','latex')
ylabel('Asymmetry parameter $\lambda$','interpreter','latex')
title('Asymmetry')
legend('GCM','1-D Modal Theory')
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);
set(gca,'linewidth',1.5)
ylim([0.5,1])
%saveas(gcf,['Plots_Proposal/Large_Domain/lambda_vs_r_',IC,'_dx',num2str(round(dx,2)),'.png'])

% Plot of magnitude of profiles

figure(3)
semilogx(R,partial_t_amp); hold on;
semilogx(R,bp_amp); hold on
semilogx(R,pb_amp); hold on;
semilogx(R,heat_amp);
xlabel('Reduction factor r','interpreter','latex')
ylabel('RMS','interpreter','latex')
title('Asymmetry')
legend('partial_t','uq_x','vq_y','heat')
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);
set(gca,'linewidth',1.5)




















% for jj = 1:length(tt)-1
%     deltat(jj) = tt(jj+1)-tt(jj);
%     
% end
% 
% for jj = 1:size(y,1)
%     
%    amp(jj) = rms(y(jj,:));
%    
% end
% 
% figure(1);
% loglog(t,amp)
% xlabel('time since last rescaling')
% ylabel('rms')
% title('Evolution of Amplitude')
% xlim([0,t(end)]);
    
