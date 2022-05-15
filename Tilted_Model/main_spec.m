% Solve the Mode-Equation as an initial value problem:
% d_tt ((rw)_xx-4*w) = (rw)_xxxx+4*w_xx
% d_tt B*w = A*w -> d_tt w = B^(-1)*A*w

clear;
close all;

% Define Grid

N = 200; 
L = 16*pi; 
x = linspace(0,L,N+1);
dx = L/(N);

k0 =  2*pi/L;
k = [0:N/2,-N/2+1:-1]*k0; k = k';

wv2 = k.^2;
wvx = k*dx;
kmax2=((N/2-1)*k0).^2;
trunc=(wv2<kmax2);
cphi = 0.715*pi;
filtr=exp(-18*(wvx-cphi).^7).*(wvx>cphi)+(wvx<=cphi);
filtr(isnan(filtr))=1;
Filtr = [filtr;filtr];

% Define time-step

t_final = 400;
tN = 16000;
dt = t_final/tN;

% define r-factor

R = [1.0;0.9;0.8;0.7;0.6;0.5;0.4;0.3;0.2;0.18;0.16;0.14;0.12;...
    0.11;0.10;0.08;0.06;0.04;0.02;0.01]; 

% define constants

beta = 0; % beta-factor
nu = 10^(-4); % hyperviscosity

% define variables

sigma = zeros(length(R),1); lambda = zeros(length(R),1); 

w_final = zeros(N+1,length(R));

[partial_t,adv_up,adv_low,adv_bp_up,adv_bp_low,adv_pb_up,adv_pb_low,heating,...
budget_w,budget_w_mode,budget_PV] = deal(zeros(length(R),1));

[w_equation,w_equation_mode,PV_equation,...
phi_equation,tau_equation] = deal(zeros(N,length(R))); 

% ode45 specifications

options = odeset('Events',@myEvent_spec);


for ii = 1:1
ii

% Initialize Values

[w_test,r_test,tau_test,phi_test] = deal(zeros(N,3));

[Phi_test] = deal(zeros(2*N,3));

% Initial Values for Phi

% sinusoidal-guess

% phi = sin(x)+sin(2*x)+sin(4*x)+sin(6*x)+sin(8*x)+sin(10*x); phi = phi(1:end-1)';
% tau = sin(x)+sin(2*x)+sin(4*x)+sin(6*x)+sin(8*x)+sin(10*x); tau = tau(1:end-1)';
% phi = fft(phi);
% tau = fft(tau);

% phi = sin(x)+sin(5*x); phi = phi(1:end-1)'; phi = fft(phi);
% tau = sin(x)+sin(5*x); tau = tau(1:end-1)'; tau = fft(tau);


% % % random-guess
%  
phi = randn(N,1); phi = phi - mean(phi); %phi = fft(phi);
tau = randn(N,1); tau = tau - mean(tau); %tau = fft(tau);

Phi0 = [tau;phi];

% Define inverse matrix

A2 = d_2x(N,dx); A2(1,:) = 1; A2inv = inv(A2);

%Propagator: 
total_time = 0;
g = [];

Phi = time_stepper_1(t_final,Phi0,R(ii),beta,nu,trunc,Filtr,N,k,dx,dt);

%[t,y] = ode45(@(t,y) Prop_beta_spectral(t,y,R(ii),beta,nu,N,k,dx),[0 t_final],Phi0);
% 
% while total_time <= t_final
% [t,y] = ode45(@(t,y) Prop_beta_spectral(t,y,R(ii),beta,nu,N,k,dx),[0 t_final],Phi0,options);
% Phi0 = y(end,:)/100;
% total_time = t(end)+total_time;
% g = cat(1,g,y);
% end

%y(end,:) = rescale(y(end,:));
dt = t(end)-t(end-1);

Phi = y(:,N+1:end)'; Tau = y(:,1:N)';

% Last step w-inversion
Phi_end = Phi(:,end); Phi_end(1) = 0;
Phi_p_end = Phi(:,end-1); Phi_p_end(1) = 0;

Tau_end = Tau(:,end); Tau_end(1) = 0;
Tau_p_end = Tau(:,end-1); Tau_p_end(1) = 0;

phi_end = A2inv*Phi_end; phi_p_end = A2inv*Phi_p_end;
tau_end = A2inv*Tau_end; tau_p_end = A2inv*Tau_p_end;

% phi = y(:,N+1:end)'; tau = y(:,1:N)';
% 
% phi_end = phi(:,end); phi_p_end = phi(:,end-1);
% tau_end = tau(:,end); tau_p_end = tau(:,end-1); 

RHS = 2*d_1x(N,dx)*Phi(:,end)+beta*d_1x(N,dx)*tau_end;
%RHS = 2*d_1x(N,dx)*Phi(:,end-1)+beta*d_1x(N,dx)*tau_p_end;

%RHS = 2*d_3x(N,dx)*phi(:,end)+beta*d_1x(N,dx)*tau(:,end);
%RHS = 2*d_3x(N,dx)*phi(:,end-1)+beta*d_1x(N,dx)*tau(:,end-1);

w = Omega_Solver(RHS,R(ii),N,dx);
r = r_factor(w,R(ii));

% final fields 

w_final(:,ii) = [w;w(1)];

% growthrate 

sigma(ii) = growthrate(y,t);
 
% % check-equation closure

w_equation(:,ii) = 1/norm(w)*(d_2x(N,dx)*(r.*w)-w-RHS);

tau_equation(:,ii) = 1/norm(Tau(:,end))*(1/dt*(Tau(:,end)-Tau(:,end-1))+d_1x(N,dx)*Phi(:,end)+beta*d_1x(N,dx)*tau_end+w);

phi_equation(:,ii) = 1/norm(Phi(:,end))*(1/dt*(Phi(:,end)-Phi(:,end-1))+d_1x(N,dx)*Tau(:,end)+beta*d_1x(N,dx)*phi_end);

% Calculate the PV-Budget

[partial_t_term,low,bp_low,pb_low,heat] = main_PV_beta(w,tau_end,tau_p_end,phi_end,phi_p_end,R(ii),beta,dt,N,dx);

%[partial_t_term,low,bp_low,pb_low,heat] = main_PV_beta_test(w,tau_end,tau_p_end,phi_end,phi_p_end,R(ii),sigma(ii),beta,dt,N,dx);

PV_equation(:,ii) = (partial_t_term + pb_low + bp_low - heat);


figure(24)
plot(x,-[heat;heat(1)]); hold on;
plot(x,[partial_t_term;partial_t_term(1)]); hold on
plot(x,[bp_low;bp_low(1)]); hold on;
plot(x,[pb_low;pb_low(1)]); hold on
plot(x,[PV_equation(:,ii);PV_equation(1,ii)])
plot(x,w_final(:,ii)); hold on;
xlabel('x'); ylabel('Magnitude')
legend('heat','Partial t','uqx','vqy','Residual','w');
title(['r=',num2str(R(ii)),' ','beta=',num2str(beta),' ','dx=',num2str(round(dx,2)),' ','sigma=',num2str(round(sigma(ii),2)),' ','rand'])
%savefig(['/net/aimsir/archive1/mkohl/Mode_Zurita_Gotor/transition/PV_Budget_rand_r',num2str(R(ii)),'_b',num2str(beta),'_dx',num2str(round(dx,2)),'.fig'])

end