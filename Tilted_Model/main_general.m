%% 2-Layer QG-Model with Latent Heating and no Internal PV-Gradients
% Tilted Interfaces (to produce zero internal PV-gradients) h_b = h_t = h = y
%and linearization around a mean shear state tau = -y produces the modelling
% equation: 
% (i) d_t phi_xx + tau_xxx - tau_x = 0
% (ii) d_t tau_xx + phi_xxx - phi_x + w = 0
% (iii) (r(w)w)_xx - w = 2*phi_xxx - phi_x
% Plot the Lower PV-Budget: PV_2 = (phi-tau)_xx + tau + h

clear;
close all;

% Define Grid

N = 200; 
L = 8*pi; 
x = linspace(0,L,N+1);
dx = L/(N);

% Define Layer Ration H1/H2

alpha = 1;

% Define time-step

t_final = 400;
tN = 16000;
dt = t_final/tN;

% define r-factor

R = [1.0,0.8,0.6,0.4,0.3,0.2,0.1,0.05,0.02,0.01];

% Initial Conditions

IC = 'rand';
%IC = 'sinusoidal';

% define variables

sigma = zeros(length(R),1); lambda = zeros(length(R),1); 

w_final = zeros(N+1,length(R));

[partial_t,adv_up,adv_low,adv_bp_up,adv_bp_low,adv_pb_up,adv_pb_low,heating,...
budget_w,budget_w_mode,budget_PV] = deal(zeros(length(R),1));

[w_equation,w_equation_mode,PV_equation,...
phi_equation,tau_equation] = deal(zeros(N,length(R))); 

% Rescale the Fields

options = odeset('Events',@myEvent);


for ii = 10:10
ii

% Initialize Values

[w_test,r_test,tau_test,phi_test] = deal(zeros(N,3));

[Phi_test] = deal(zeros(2*N,3));

% Initial Conditions: Choose between random and sinusoidal guess

if strcmp(IC,'rand') == 1
% random-guess
phi =   randn(1,N); phi = phi - mean(phi); 
tau =   randn(1,N); tau = tau - mean(tau); 
else

% sinusoidal-guess

phi = sin(x)+sin(2*x)+sin(4*x)+sin(6*x)+sin(8*x)+sin(10*x); phi = phi(1:end-1);
tau = sin(x)+sin(2*x)+sin(4*x)+sin(6*x)+sin(8*x)+sin(10*x); tau = tau(1:end-1);

%phi = sin(5*x); phi = phi(1:end-1);
%tau = sin(5*x); tau = tau(1:end-1);

end

Phi0 = [tau,phi]';

% Define inverse matrix

A2 = d_2x(N,dx); A2(1,:) = 1; A2inv = inv(A2);

%Propagator: 
total_time = 0;
g = [];

while total_time <= t_final
[t,y] = ode45(@(t,y) Prop(t,y,A2inv,R(ii),N,dx),[0 t_final],Phi0,options);
Phi0 = y(end,:)/100;
total_time = t(end)+total_time;
g = cat(1,g,y);
end

dt = t(end)-t(end-1);

Phi = y(:,N+1:end)'; Tau = y(:,1:N)';

% Last step w-inversion
Phi_end = Phi(:,end); Phi_end(1) = 0;
Phi_p_end = Phi(:,end-1); Phi_p_end(1) = 0;

Tau_end = Tau(:,end); Tau_end(1) = 0;
Tau_p_end = Tau(:,end-1); Tau_p_end(1) = 0;

phi_end = A2inv*Phi_end; phi_p_end = A2inv*Phi_p_end;
tau_end = A2inv*Tau_end; tau_p_end = A2inv*Tau_p_end;

RHS = 2*d_1x(N,dx)*Phi(:,end)-d_1x(N,dx)*phi_end;
RHS1 = 2*d_1x(N,dx)*Phi(:,end);
RHS2 = -d_1x(N,dx)*phi_end;

w = Omega_Solver(RHS,R(ii),N,dx);
r = r_factor(w,R(ii));

% final fields 

w_final(:,ii) = [w;w(1)];

% growthrate 

sigma(ii) = growthrate(y,t);

% Calculate Asymmetry

lambda(ii) = asymmetry(w_final(:,ii));
 
% % check-equation closure

w_equation(:,ii) = 1/norm(w)*(d_2x(N,dx)*(r.*w)-w-RHS);

tau_equation(:,ii) = 1/norm(Tau(:,end))*(1/dt*(Tau(:,end)-Tau(:,end-1))+d_1x(N,dx)*Phi(:,end)-d_1x(N,dx)*phi_end+w);

phi_equation(:,ii) = 1/norm(Phi(:,end))*(1/dt*(Phi(:,end)-Phi(:,end-1))+d_1x(N,dx)*Tau(:,end)-d_1x(N,dx)*tau_end);

% Calculate the PV-Budget

%[partial_t_term,low,bp_low,pb_low,heat] = main_PV_Budget(w,tau_end,tau_p_end,phi_end,phi_p_end,R(ii),dt,N,dx);

[partial_t_term,low,bp_low,pb_low,heat] = main_PV_Budget_sigma(w,tau_end,phi_end,R(ii),sigma(ii),N,dx);

[partial_t_term_up,up,bp_up,pb_up,heat_up] = main_PV_Budget_sigma(w,tau_end,phi_end,R(ii),sigma(ii),N,dx);

PV_equation(:,ii) = (partial_t_term + pb_low + bp_low - heat);


% figure(24)
% plot(x,-[heat;heat(1)]); hold on;
% plot(x,[partial_t_term;partial_t_term(1)]); hold on
% plot(x,[bp_low;bp_low(1)]); hold on;
% plot(x,[pb_low;pb_low(1)]); hold on
% plot(x,[PV_equation(:,ii);PV_equation(1,ii)])
% plot(x,w_final(:,ii)); hold on;
% xlabel('x'); ylabel('Magnitude')
% ax = gca;
% ax.XTick = [0,2*pi,4*pi,6*pi,8*pi];
% ax.XTickLabel = {'0','2\pi','4\pi','6\pi','8\pi'};
% title('w-profile for r=0.1')
% set(gca,'fontsize', 14);
% xlim([0 8*pi])
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% set(gca,'linewidth',1.5)
% legend('heat','Partial t','uqx','vqy','Residual','w');
% %title(['r=',num2str(R(ii)),' ','beta=',num2str(beta),' ','dx=',num2str(round(dx,2)),' ','sigma=',num2str(round(sigma(ii),2)),' ','rand'])
% title(['r=',num2str(R(ii))])
%saveas(gcf,['Plots_Proposal_Main/Large_Domain/PV_Budget_',IC,'_r',num2str(R(ii)),'_b',num2str(beta),'_dx',num2str(round(dx,2)),'.png'])
%close all;

%savefig(['/net/aimsir/archive1/mkohl/Mode_Zurita_Gotor/transition/PV_Budget_rand_r',num2str(R(ii)),'_b',num2str(beta),'_dx',num2str(round(dx,2)),'.fig'])

%qx = 2/sigma(ii)*d_1x(N,dx)*partial_t_term;
%qx = 2*partial_t_term;
%omega = d_2x(N,dx)*(r.*w)-w;

% LHS = r.*d_2x(N,dx)*w-w;
% LHS1 = r.*d_2x(N,dx)*w;
% LHS2 = -w;
% 
% figure(5)
% plot(LHS); hold on;
% plot(LHS1); hold on;
% plot(LHS2)
% 
% figure(20)
% plot(x,[RHS;RHS(1)]); hold on;
% plot(x,[RHS1;RHS1(1)]); hold on;
% plot(x,[RHS2;RHS2(1)]); hold on;
% plot(x,[qx;qx(1)]);
% %plot(x,[omega;omega(1)]);
% plot(x,[w;w(1)]);
% legend('RHS','2*phi_{xxx}','-phi_{x}','2*qx','w')
% 
% figure(24)
% plot(-[heat;heat(1)]); hold on;
% plot([partial_t_term;partial_t_term(1)]); hold on
% plot([bp_low;bp_low(1)]); hold on;
% plot([pb_low;pb_low(1)]); hold on
% plot([PV_equation(:,ii);PV_equation(1,ii)])
% plot(w_final(:,ii)); hold on;
% ylabel('Magnitude')
% title('w-profile for r=0.1')
% set(gca,'fontsize', 14);
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% set(gca,'linewidth',1.5)
% legend('heat','Partial t','uqx','vqy','Residual','w');
% %title(['r=',num2str(R(ii)),' ','beta=',num2str(beta),' ','dx=',num2str(round(dx,2)),' ','sigma=',num2str(round(sigma(ii),2)),' ','rand'])
% title(['r=',num2str(R(ii))])

partial_t_term = [partial_t_term;partial_t_term(1)]/max(bp_low);
pb_low = [pb_low;pb_low(1)]/max(bp_low);
heat = [heat;heat(1)]/max(bp_low);
Residual = [PV_equation(:,ii);PV_equation(1,ii)]/max(bp_low);
bp_low = [bp_low;bp_low(1)]/max(bp_low);

index = find(max(bp_low)==bp_low);
up = index+60;
low = index-60;

figure(24)
plot(x(low:up),partial_t_term(low:up),'linewidth',1.6); hold on
plot(x(low:up),bp_low(low:up),'linewidth',1.6); hold on;
plot(x(low:up),pb_low(low:up),'linewidth',1.6); hold on
plot(x(low:up),-heat(low:up),'linewidth',1.6); hold on;
plot(x(low:up),Residual(low:up),'linewidth',1.6)
%plot(x,w_final(:,ii),'linewidth',1.6); hold on;
xlabel('x'); ylabel('Magnitude')
xlim([x(low) x(up)])
title('w-profile for r=0.1')
set(gca,'fontsize', 12);
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'linewidth',1.5)
legend({'\sigma q','uq_{x}','vq_{y}','-LH','Sum'},'Location','south','Orientation','horizontal');
legend boxoff
title('PV-Budget')
%saveas(gcf,'/disk7/mkohl/Mode_Zurita_Gotor/Plots_Proposal_Main/Large_Domain/PV_Budget_DRV_r001_HD','epsc')
close all;
% % 
figure(30)
plot(x,w_final(:,ii)/norm(w_final(:,ii)),'linewidth',1.6); hold on;
xlabel('x'); ylabel('w')
% ax = gca;
% ax.XTick = [0,2*pi,4*pi,6*pi,8*pi];
% ax.XTickLabel = {'0','2\pi','4\pi','6\pi','8\pi'};
title('r=0.01')
set(gca,'fontsize', 12);
xlim([0 L])
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'linewidth',1.5)
%saveas(gcf,'/disk7/mkohl/Mode_Zurita_Gotor/Plots_Proposal_Main/Large_Domain/w_DRV_beta0_r001_HD','epsc')
end