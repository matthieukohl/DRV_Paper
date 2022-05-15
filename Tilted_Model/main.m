%% Tilted 2-Layer QG-Model with Latent Heating and no Internal PV-Gradients
% Tilted Interfaces (to produce zero internal PV-gradients) h_b = h_t = h = y
% and linearization around a mean shear state tau = -y produces the modelling
% equation: 
% (i) d_t phi_xx + tau_xxx - h*tau_x = 0
% (ii) d_t tau_xx + phi_xxx - h*phi_x + w = 0
% (iii) (r(w)w)_xx - w = 2*phi_xxx - h*phi_x
% h = 0 for the classic untilted model; h = 1 for the tilted model

clear;
close all;


% Define Grid

%N = 2000;
%L = 16*pi;
N = 200;
L = 8*pi;
x = linspace(0,L,N+1);
dx = L/(N);


% Define time-step

t_final = 200;
% tN = 16000;
% dt = t_final/tN;

% define r-factor

R = [1,0.9,0.8,0.6,0.4,0.3,0.2,0.1,0.05,0.02,0.01];

%R = [1,0.8,0.6,0.4,0.35,0.3,0.25,0.22,0.21,0.2,0.19,0.18,0.15,0.1,0.08,0.06,0.04,0.01,0.005,0.001];

%R = [0.9,0.8,0.6,0.4,0.2,0.1,0.05,0.01,0.005,0.001];

%R = 0.5;


% untilted h = 0, tilted h = 1;

h = 0;

% Initial Conditions

IC = 'rand';
%IC = 'sinusoidal';

% define variables

sigma = zeros(length(R),1); lambda = zeros(length(R),1); 

b = zeros(length(R),1);

w_final = zeros(N+1,length(R));

[partial_t,adv_up,adv_low,adv_bp_up,adv_bp_low,adv_pb_up,adv_pb_low,heating,...
budget_w,budget_w_mode,budget_PV] = deal(zeros(length(R),1));

[w_equation,w_equation_mode,PV_equation,...
phi_equation,tau_equation] = deal(zeros(N,length(R))); 

% Rescale the Fields

options = odeset('Events',@myEvent);


for ii = 1:length(R)
ii
tic

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

%phi = sin(x)+sin(2*x)+sin(4*x)+sin(6*x)+sin(8*x)+sin(10*x); phi = phi(1:end-1);
%tau = sin(x)+sin(2*x)+sin(4*x)+sin(6*x)+sin(8*x)+sin(10*x); tau = tau(1:end-1);

phi = sin(x); phi = phi(1:end-1);
tau = sin(x); tau = tau(1:end-1);

end

Phi0 = [tau,phi]';

% Define inverse matrix

A2 = d_2x(N,dx); A2(1,:) = 1; A2inv = inv(A2);

%Propagator: 
total_time = 0;
g = [];
W = [];
s_mode = [];
time = 0;
nrescale = 0;

while total_time < t_final
[t,y] = ode45(@(t,y) Prop(t,y,A2inv,R(ii),h,N,dx),[0 t_final-total_time],Phi0,options);
Phi0 = y(end,:)/100;
total_time = t(end)+total_time;
g = cat(1,g,y);

% calculate the growthrate over the time of a rescaling event
% s = 1/(t(end)-t(1))*log(rms(y(end,:))/rms(y(1,:)));
% 
% s_mode = cat(1,s_mode,s);

% calculate the growthrate for every timestep

time = cat(1,time,t(2:end)+time(end));

s = zeros(length(t),1);
for kk = 2:length(t)
    s(kk) = 1/(t(kk)-t(kk-1))*log(rms(y(kk,:))/rms(y(kk-1,:)));

end
s_mode = cat(1,s_mode,s(2:end));


nrescale = nrescale +1;
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

RHS = 2*d_1x(N,dx)*Phi(:,end)-h*d_1x(N,dx)*phi_end;
RHS1 = 2*d_1x(N,dx)*Phi(:,end);
RHS2 = -d_1x(N,dx)*phi_end;

w = Omega_Solver(RHS,R(ii),N,dx);
r = r_factor(w,R(ii));

% final fields 

w_final(:,ii) = [w;w(1)];

% growthrate: average the growth rate over the last five rescale events
% sigma(ii) = mean(s_mode(end-4:end));

% growthrate: average the stepwise growthrate over the last 5 simulation
% steps
[val,index] = min(abs(time-195));
sigma(ii) = mean(s_mode(index:end));

% Calculate Asymmetry

lambda(ii) = asymmetry(w_final(:,ii));

% Calculate Half-Ascent Area

index = find(w>0); 

b(ii) = (length(index)+1)*dx/2;


% % check-equation closure

% w_equation(:,ii) = 1/norm(w)*(d_2x(N,dx)*(r.*w)-w-RHS);
% 
% tau_equation(:,ii) = 1/norm(Tau(:,end))*(1/dt*(Tau(:,end)-Tau(:,end-1))+d_1x(N,dx)*Phi(:,end)-d_1x(N,dx)*phi_end+w);
% 
% phi_equation(:,ii) = 1/norm(Phi(:,end))*(1/dt*(Phi(:,end)-Phi(:,end-1))+d_1x(N,dx)*Tau(:,end)-d_1x(N,dx)*tau_end);

% Calculate the PV-Budget

%[partial_t_term,low,bp_low,pb_low,heat] = main_PV_Budget(w,tau_end,tau_p_end,phi_end,phi_p_end,R(ii),dt,N,dx);

[partial_t_term,low,bp_low,pb_low,heat] = main_PV_Budget_sigma(w,tau_end,phi_end,R(ii),sigma(ii),N,dx);

[partial_t_term_up,up,bp_up,pb_up,heat_up] = main_PV_Budget_up_sigma(w,tau_end,phi_end,R(ii),sigma(ii),N,dx);

pb_low = (-1+h)*d_1x(N,dx)*(phi_end-tau_end);

PV_equation(:,ii) = (partial_t_term + pb_low + bp_low - heat);

PV_equation_up = (partial_t_term_up + pb_up + bp_up - heat_up);


partial_t_term = [partial_t_term;partial_t_term(1)]/max(bp_low);
pb_low = [pb_low;pb_low(1)]/max(bp_low);
heat = [heat;heat(1)]/max(bp_low);
Residual = [PV_equation(:,ii);PV_equation(1,ii)]/max(bp_low);
bp_low = [bp_low;bp_low(1)]/max(bp_low);

partial_t_term_up = [partial_t_term_up;partial_t_term_up(1)]/(abs(min(bp_up)));
pb_up = [pb_up;pb_up(1)]/(abs(min(bp_up)));
heat_up = [heat_up;heat_up(1)]/(abs(min(bp_up)));
Residual_up = [PV_equation_up;PV_equation_up(1)]/(abs(min(bp_up)));
bp_up = [bp_up;bp_up(1)]/(abs(min(bp_up)));

index = find(max(bp_low)==bp_low);
up = index+60;
low = index-30;

% % figure(24)
%  plot(x(low:up),partial_t_term(low:up),'linewidth',1.6); hold on
%  plot(x(low:up),-bp_low(low:up),'linewidth',1.6); hold on;
%  plot(x(low:up),-pb_low(low:up),'linewidth',1.6); hold on
%  plot(x(low:up),heat(low:up),'linewidth',1.6); hold on;
%  plot(x(low:up),Residual(low:up),'linewidth',1.6)
% % %plot(x,w_final(:,ii),'linewidth',1.6); hold on;
%  xlabel('x'); ylabel('Magnitude')
%  xlim([x(low) x(up)])
%  ylim([-1.5 1.5])
%  title('w-profile for r=0.1')
%  set(gca,'fontsize', 12);
%  set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
%  set(gca,'linewidth',1.5)
%  legend({'\sigma q','uq_{x}','vq_{y}','-LH','Sum'},'Location','south','Orientation','horizontal');
%  legend boxoff
%  title('PV-Budget Low')
%saveas(gcf,'/disk7/mkohl/Mode_Kohl_OGorman/figures_work_in_progress/PV_Budget_DRV_q2_r001_HD','epsc')
%close all;


% figure(25)
% plot(x(low:up),partial_t_term_up(low:up),'linewidth',1.6); hold on
% plot(x(low:up),bp_up(low:up),'linewidth',1.6); hold on;
% plot(x(low:up),pb_up(low:up),'linewidth',1.6); hold on
% plot(x(low:up),-heat_up(low:up),'linewidth',1.6); hold on;
% plot(x(low:up),Residual_up(low:up),'linewidth',1.6)
% %plot(x,w_final(:,ii),'linewidth',1.6); hold on;
% xlabel('x'); ylabel('Magnitude')
% xlim([x(low) x(up)])
% set(gca,'fontsize', 12);
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% set(gca,'linewidth',1.5)
% legend({'\sigma q','uq_{x}','vq_{y}','-LH','Sum'},'Location','north','Orientation','horizontal');
% legend boxoff
% title('PV-Budget Up')
%saveas(gcf,'/disk7/mkohl/Mode_Zurita_Gotor/Plots_Proposal_Main/Large_Domain/PV_Budget_DRV_q1_r001_HD','epsc')
%close all;


% figure(30)
% plot(x,w_final(:,ii)/norm(w_final(:,ii)),'linewidth',1.6); hold on;
% xlabel('x'); ylabel('w')
% % ax = gca;
% % ax.XTick = [0,2*pi,4*pi,6*pi,8*pi];
% % ax.XTickLabel = {'0','2\pi','4\pi','6\pi','8\pi'};
% title(['r=',num2str(R(ii))])
% set(gca,'fontsize', 12);
% xlim([0 L])
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% set(gca,'linewidth',1.5)
%saveas(gcf,'/disk7/mkohl/Mode_Zurita_Gotor/Plots_Proposal_Main/Large_Domain/w_DRV_beta0_r001_q1q2_HD','epsc')
toc
end

%save('DRV_Analytical/Time_Marching_L16pi_N2000','R','N','L','sigma','b')

%save('/disk7/mkohl/Asymmetry_Reanalysis2/data_theory/modal_theory_sh_N800_L32pi.mat','lambda');

%save sample DRV
%R = R(ii); sigma = sigma(ii);
%save('DRV_Paper/DRV_example_L16pi_N2000.mat','Phi','Tau','phi_end','tau_end','w','R','sigma','x','L','N','dx');


