% Modal Instability calculations with tilted boundaries at both
% the top h1 and bottom h2

clear;
close all;


% Define Grid

N = 200;
L = 8*pi; 
x = linspace(0,L,N+1);
dx = L/(N);


% Define time-step

%t_final = 600;
t_final = 200;
tN = 16000;
dt = t_final/tN;

% define r-factor


%R = [0.01,0.05,0.1,0.3,0.5,0.7,1.0];

R = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];



% Initial Conditions

IC = 'rand';
%IC = 'sinusoidal';

% beta-factor

beta = zeros(size(R));

% alpha - drag

alpha1 = zeros(length(R));
alpha2 = zeros(length(R));
alpha3 = zeros(length(R)); 

% tilted boundaries h1 and h2

%h = 0:0.1:2; 

h = [0:0.1:0.3,0.35,0.4:0.1:2];

h1 = h; h2 = h;

% define variables

[sigma,sigma2,b,mode_or_drv,mode_or_drv2] = deal(zeros(length(R),length(h))); 

ww = zeros(length(R),length(h),N);
phi_stream = zeros(length(R),length(h),N);
tau_stream = zeros(length(R),length(h),N);

w_final = zeros(N+1,length(R));

[partial_t,adv_up,adv_low,adv_bp_up,adv_bp_low,adv_pb_up,adv_pb_low,heating,...
budget_w,budget_w_mode,budget_PV] = deal(zeros(length(R),1));

[w_equation,w_equation_mode,PV_equation,...
phi_equation,tau_equation] = deal(zeros(N,length(R))); 

% Rescale the Fields

options = odeset('Events',@myEvent);


for ii = 1:length(R) %length(R) %length(R)
  
 for jj = 5:length(h) %length(h) %length(h):length(h) 
     tic
ii
jj

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
s_mode = [];
time = 0;
g = [];
W = [];

tt = [];
yy = [];

while total_time < t_final
[t,y] = ode45(@(t,y) Prop(t,y,A2inv,R(ii),N,dx,h(jj),h(jj),beta(ii),alpha1(ii),alpha2(ii),alpha3(ii)),[0 t_final-total_time],Phi0,options);

Phi0 = y(end,:)/100;
total_time = t(end)+total_time;
g = cat(1,g,y);

% % calculat the growthrate over the time of a rescaling event
% s = 1/(t(end)-t(1))*log(rms(y(end,:))/rms(y(1,:)));
% 
% s_mode = cat(1,s_mode,s);

time = cat(1,time,t(2:end)+time(end));

s = zeros(length(t),1);
for kk = 2:length(t)
    s(kk) = 1/(t(kk)-t(kk-1))*log(rms(y(kk,:))/rms(y(kk-1,:)));

end
s_mode = cat(1,s_mode,s(2:end));

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


RHS = 2*d_1x(N,dx)*Phi(:,end)-1/2*(h1(jj)+h2(jj))*d_1x(N,dx)*phi_end + 1/2*(h2(jj)-h1(jj))*d_1x(N,dx)*tau_end...
+beta(ii)*d_1x(N,dx)*tau_end-(alpha3(ii)-alpha2(ii))*d_2x(N,dx)*tau_end;


RHS1 = 2*d_1x(N,dx)*Phi(:,end);
RHS2 = -d_1x(N,dx)*phi_end;

w = Omega_Solver(RHS,R(ii),N,dx);
r = r_factor(w,R(ii));

% w-field

ww(ii,jj,:) = w;

% streamfunctions

phi_stream(ii,jj,:) = phi_end;
tau_stream(ii,jj,:) = tau_end;

% growthrate 

% sigma(ii,jj) = growthrate(y,t);
% 
% sigma2(ii,jj) = 1/(t(end)-t(1))*log(rms(y(end,:))/(rms(y(1,:))));

% sigma(ii,jj) = mean(s_mode(end-4:end));

[val,index] = min(abs(time-(t_final-5)));
sigma(ii,jj) = mean(s_mode(index:end));

% Calculate Half-Ascent Area

b(ii,jj) = dx*(min(abs(find(w==max(w))-find(w<0))));

% Mode or DRV - that is the question; DRV if single local maximum

ind = islocalmax(w); 
ind(w<0) = 0; % sometimes DRV solutions have a local max when w<0


if sum(ind)==1
 mode_or_drv(ii,jj) = 1;
else
 mode_or_drv(ii,jj) = -1;
end

% plot(x(1:end-1),w); xlim([0 x(end-1)]); hold on;
% title(['\sigma=',num2str(sigma(ii,jj))]);


% Calculate the PV-Budget

%[partial_t_term,low,bp_low,pb_low,heat] = main_PV_Budget(w,tau_end,tau_p_end,phi_end,phi_p_end,R(ii),dt,N,dx);

[partial_t_term,low,bp_low,pb_low,heat] = main_PV_Budget_low_sigma(w,tau_end,phi_end,R(ii),sigma(ii,jj),N,dx,h1(jj),h2(jj),beta(ii));

% [partial_t_term,low,bp_low,pb_low,heat] = main_PV_Budget(w,tau_end,tau_p_end,phi_end,phi_p_end,R(ii),dt,N,dx,beta(ii));
% 
[partial_t_term_up,up,bp_up,pb_up,heat_up] = main_PV_Budget_up_sigma(w,tau_end,phi_end,R(ii),sigma(ii,jj),N,dx,h1(jj),h2(jj),beta(ii));

% 

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
up = index+30;
low = index-30;

% figure(24)
 plot(x(low:up),partial_t_term(low:up),'linewidth',1.6); hold on
 plot(x(low:up),-bp_low(low:up),'linewidth',1.6); hold on;
 plot(x(low:up),-pb_low(low:up),'linewidth',1.6); hold on
 plot(x(low:up),heat(low:up),'linewidth',1.6); hold on;
 plot(x(low:up),Residual(low:up),'linewidth',1.6)
% %plot(x,w_final(:,ii),'linewidth',1.6); hold on;
 xlabel('x'); ylabel('Magnitude')
 xlim([x(low) x(up)])
 ylim([-1.5 1.5])
 title('w-profile for r=0.1')
 set(gca,'fontsize', 12);
 set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
 set(gca,'linewidth',1.5)
 legend({'\sigma q','uq_{x}','vq_{y}','-LH','Sum'},'Location','south','Orientation','horizontal');
 legend boxoff
 title('PV-Budget Low')
%saveas(gcf,'/disk7/mkohl/Mode_Kohl_OGorman/figures_work_in_progress/PV_Budget_DRV_q2_r001_HD','epsc')
%close all;




% figure(1)
% plot(x,[partial_t_term_up;partial_t_term_up(1)]+2*max(abs(partial_t_term_up)),'Linewidth',2.5,'Color','b'); hold on;
% plot(x,[partial_t_term; partial_t_term(1)],'Linewidth',2.5,'Color','r'); 
% set(gca,'fontsize', 14);
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% title('PV Anomalies')
% legend('Up','Low','Location','East'); legend boxoff


% figure(24)
% plot(-[heat;heat(1)]); hold on;
% plot([partial_t_term;partial_t_term(1)]); hold on
% plot([bp_low;bp_low(1)]); hold on;
% plot([pb_low;pb_low(1)]); hold on
% plot([PV_equation(:,ii);PV_equation(1,ii)])
% %plot(w_final(:,ii)); hold on;
% ylabel('Magnitude')
% set(gca,'fontsize', 14);
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% set(gca,'linewidth',1.5)
% legend('heat','Partial t','uqx','vqy','Residual','w');
% %legend('heat','Partial t','vqy','Residual');
% legend('boxoff')
% %title(['r=',num2str(R(ii)),' ','beta=',num2str(beta),' ','dx=',num2str(round(dx,2)),' ','sigma=',num2str(round(sigma(ii),2)),' ','rand'])
% title(['r=',num2str(R(ii)),'','\beta=',num2str(beta(ii))])

% figure(25)
% plot(x,w_final(:,ii),'Linewidth',2.5); hold on;
% ylabel('Magnitude')
% set(gca,'fontsize', 14);
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% set(gca,'linewidth',1.5)
% title(['r=',num2str(round(R(ii),2)),' ','\beta=',num2str(round(beta(ii),2)),...
%     '  ','\sigma=',num2str(round(sigma(ii),2))])
% legend('heat','Partial t','uqx','vqy','Residual','w');
% legend('boxoff')
% title(['r=',num2str(R(ii)),' ','beta=',num2str(beta),' ','dx=',num2str(round(dx,2)),' ','sigma=',num2str(round(sigma(ii),2)),' ','rand'])
%  
toc
end

end

%save('/disk7/mkohl/Mode_Kohl_OGorman_beta_no_gradients/mode_or_drv/mode_or_drv_N200_L8pi_t200_with_stream.mat','L','N','t_final','R','h','ww','phi_stream','tau_stream','sigma','sigma2',...
 %  'b','mode_or_drv')

