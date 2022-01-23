%% 2-Layer QG-Model with Latent Heating and no Internal PV-Gradients
% Tilted Interfaces (to produce zero internal PV-gradients) h_b = h_t = h = y
% and linearization around a mean shear state tau = -y produces the modelling
% equation: 
% (i) d_t phi_xx + tau_xxx - tau_x = 0
% (ii) d_t tau_xx + phi_xxx - phi_x + w = 0
% (iii) (r(w)w)_xx - w = 2*phi_xxx - phi_x
% Plot the Lower PV-Budget: PV_2 = (phi-tau)_xx + tau + h

clear;
close all;


% Define Grid

N = 400;
L = 16*pi; 
x = linspace(0,L,N+1);
dx = L/(N);


% Define time-step

t_final = 200;
tN = 16000;
dt = t_final/tN;

% define r-factor

%load('GCM_Modal.mat')

load('2D_3D_data.mat')

%R = R_eff; R = R_mean;

R = R_mid;

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

% Calculate Half-Ascent Area

index = find(w>0); 

b(ii) = (length(index)+1)*dx/2;


% % check-equation closure

w_equation(:,ii) = 1/norm(w)*(d_2x(N,dx)*(r.*w)-w-RHS);

tau_equation(:,ii) = 1/norm(Tau(:,end))*(1/dt*(Tau(:,end)-Tau(:,end-1))+d_1x(N,dx)*Phi(:,end)-d_1x(N,dx)*phi_end+w);

phi_equation(:,ii) = 1/norm(Phi(:,end))*(1/dt*(Phi(:,end)-Phi(:,end-1))+d_1x(N,dx)*Tau(:,end)-d_1x(N,dx)*tau_end);

% Calculate the PV-Budget

%[partial_t_term,low,bp_low,pb_low,heat] = main_PV_Budget(w,tau_end,tau_p_end,phi_end,phi_p_end,R(ii),dt,N,dx);

[partial_t_term,low,bp_low,pb_low,heat] = main_PV_Budget_sigma(w,tau_end,phi_end,R(ii),sigma(ii),N,dx);

[partial_t_term_up,up,bp_up,pb_up,heat_up] = main_PV_Budget_up_sigma(w,tau_end,phi_end,R(ii),sigma(ii),N,dx);


PV_equation(:,ii) = (partial_t_term + pb_low + bp_low - heat);

PV_equation_up = (partial_t_term_up + pb_up + bp_up - heat_up);


end

sigma_DRV_mid = sigma;
b_DRV_mid = b;
save('/disk7/mkohl/Mode_Kohl_OGorman/2D_3D_comparison/DRV_mid.mat','sigma_DRV_mid',...
'b_DRV_mid')

