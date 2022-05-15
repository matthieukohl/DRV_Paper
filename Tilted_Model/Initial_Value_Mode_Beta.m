% Solve the Mode-Equation with beta as an initial value problem:
% d_tt ((rw)_xxxx-4*w_xx)+d_t(2*alpha*(rw)_xxx-4*alpha*w_x) =
% (rw)_xxxxxx+4*w_xxxx-alpha*(rw)_xx
% d_tt A(w)*w + d_t B(w)*w = C(w)*w

clear;
close all;

% Define Grid

N = 200; 
L = 8*pi; 
x = linspace(0,L,N+1);
dx = L/(N);

% Define time-step

t_final = 200;
tN = 16000;
dt = t_final/tN;

% Define constants

alpha = 1; % alpha = beta*L_D^2/U

% define Wave-numbers for further analysis

k0 =  2*pi/L;
k = [0:N/2,-N/2+1:-1]*k0;

% define r-factor

R = [1.0;0.9;0.8;0.7;0.6;0.5;0.4;0.3;0.2;0.18;0.16;0.14;0.12;...
    0.11;0.10;0.08;0.06;0.04;0.02;0.01]; 

% define variables

sigma = zeros(length(R),1); lambda = zeros(length(R),1); 

w_final = zeros(N+1,length(R));

[partial_t,adv_up,adv_low,adv_bp_up,adv_bp_low,adv_pb_up,adv_pb_low,heating,...
budget_w,budget_w_mode,budget_PV] = deal(zeros(length(R),1));

[w_equation,w_equation_mode,equation_one,equation_two,PV_equation,...
    phi_equation,tau_equation] = deal(zeros(N,length(R))); 

% Initialize Solution

w_new = sin(x)+sin(2*x)+sin(4*x)+sin(6*x)+sin(8*x)+sin(10*x); 
w_new = w_new(1:end-1);


for ii = 20:20
ii

% Initialize Values

[w_test,u_test,r_test,phi_test,tau_test] = deal(zeros(N,3));

% Initial values    
w = randn(1,N); %sin(x)+sin(2*x)+sin(4*x)+sin(6*x)+sin(8*x)+sin(10*x); w = w(1:end-1);
w_t = zeros(1,N);  
w = [w,w_t]';  

%Initial Values for Phi

phi = zeros(1,N);
tau = zeros(1,N);
Phi = [tau,phi]';

%Propagator for tau, phi fields
A1 = d_1x_4(N,dx);
A2 = d_2x_4(N,dx); A2(1,:) = 1; A2_inv = inv(A2);
A3 = d_3x_4(N,dx);

M = -[alpha*A2_inv*A1,A2_inv*A3;A2_inv*A3,alpha*A2_inv*A1];

Propagator = diag(ones(1,2*N))-dt*M; 

count = 0;
index = 0;
time = 0;
    
for tt = 1:tN
tt
count = count +1;   
% Define r-factor

r = r_factor(w(1:N),R(ii));

% Define Propagator C = [0,1;0,B^(-1)*A]

A = A_beta(r,N,dx); A(1,:) = 1;
B = B_beta(alpha,r,N,dx); B(1,:) = 1;
C = C_beta(alpha,r,N,dx); C(1,:) = 0;
 
% alternative time differencing

M = [-inv(A)*B,inv(A);C,zeros(N,N)];

D = diag(ones(1,2*N))-dt*M; 

% Propagate solution forward using implicit euler 

w = D\w;

av(tt) = mean(w(1:N)/norm(w(1:N)));

% Propagate tau and phi fields forward using implicit euler


c = [-2*w(1:N);zeros(N,1)];

Phi = inv(Propagator)*(Phi+dt*c);


% %Plot evolution of w
% if mod(tt,200)==0
% plot(x(1:N),w(1:N)); 
% pause(0.1)
% hold off;
% mean(w(1:N))
% end

% Store last three w-values to check validity of equation

if tt == tN-2 || tt == tN-1 || tt == tN
   
   index = index+1;
   w_test(:,index) = w(1:N);
   u_test(:,index)= w(N+1:end);
   r_test(:,index) = r_factor(w(1:N),R(ii));
   tau_test(:,index) = Phi(1:N);
   phi_test(:,index) = Phi(N+1:end);
   Phi_test(:,index) = Phi;
   
       
end

end

% Store final w-fields and the asymmetry factor

w_final(:,ii) = [w(1:N);w(1)];
lambda(ii) = asymmetry(w_final(:,ii)/norm(w_final(:,ii)));

% Calculate growthrate assuming exponential increase at late time

sigma(ii) = 1/dt*log(mean(w_test(:,3)./w_test(:,2)));

% Calculate residual of d_tt ((rw)_xx-w) = ((rw)_4x+w_2x)
% using backward Euler

RHS  = C_beta(alpha,r_test(:,3),N,dx)*w_test(:,3);

LHS_tt = A_beta(r_test(:,3),N,dx)*w_test(:,3);
LHSm1_tt = A_beta(r_test(:,2),N,dx)*w_test(:,2);
LHSm2_tt = A_beta(r_test(:,1),N,dx)*w_test(:,1);

LHS_t = B_beta(alpha,r_test(:,3),N,dx)*w_test(:,3);
LHSm1_t = B_beta(alpha,r_test(:,2),N,dx)*w_test(:,2);

% How well does the w-equation close?

equation_one(:,ii) = 1/norm(w_test(:,3))*(w(N+1:end)-B_beta(alpha,r_test(:,3),N,dx)*w_test(:,3)...
    -(A_beta(r_test(:,3),N,dx)*w_test(:,3)-A_beta(r_test(:,2),N,dx)*w_test(:,2))/dt);

equation_two(:,ii) = 1/norm(w_test(:,3))*((u_test(:,3)-u_test(:,2))/dt-C_beta(alpha,r_test(:,3),N,dx)*w_test(:,3));

w_equation(:,ii) = 1/(norm(w_test(:,3)))*((LHS_tt+LHSm2_tt-2*LHSm1_tt)/dt^2 + (LHS_t-LHSm1_t)/dt - RHS);

% Calculate tau and phi

[PV_tendency,low,bp_low,pb_low,beta_low,heat] = PV_Budget_Beta(w_test(:,3),tau_test(:,3),tau_test(:,2),phi_test(:,3),phi_test(:,2),R(ii),alpha,dt,N,dx);

PV_equation(:,ii) = PV_tendency+bp_low+pb_low+beta_low-heat;

%Plot w-solution

figure(22)
plot(x,w_final(:,ii)); hold on;
xlabel('x'); ylabel('Magnitude')
legend('w');
title(['r=',num2str(R(ii)),' ','and \beta=',num2str(alpha),' ','and \sigma=',num2str(sigma(ii)/2)])
%xlim([0 4*pi])
%Plot w-residual

figure(23)
plot(x,[w_equation(:,ii);w_equation(1,ii)]); hold on;
xlabel('x'); ylabel('Magnitude')
legend('w-Residual');
xlim([0 4*pi])

% Plot different terms

figure(24)
plot(x,-[heat;heat(1)]); hold on;
plot(x,[PV_tendency;PV_tendency(1)]); hold on
plot(x,[bp_low;bp_low(1)]); hold on;
plot(x,[pb_low;pb_low(1)]); hold on
plot(x,[beta_low;beta_low(1)]); hold on
plot(x,[PV_equation(:,ii);PV_equation(1,ii)])
plot(x,w_final(:,ii)); hold on;
xlabel('x'); ylabel('Magnitude')
legend('heat','Partial t','uqx','vqy','beta','Residual','w');
title(['PV-Budget Lower Level r=',num2str(R(ii)),' ','and beta=',num2str(alpha)])
xlim([0 4*pi])

% Spectral

% figure(24)
% plot(x,-[heat;heat(1)]); hold on;
% plot(x,[PV_tendency,PV_tendency(1)]); hold on
% plot(x,[adv_bp,adv_bp(1)]); hold on;
% plot(x,[adv_pb,adv_pb(1)]); hold on
% plot(x,[adv_beta,adv_beta(1)]); hold on
% plot(x,[PV_equation(:,ii);PV_equation(1,ii)])
% plot(x,w_final(:,ii)); hold on;
% xlabel('x'); ylabel('Magnitude')
% legend('heat','Tendency','uqx','vqy','beta','Residual','w');
% title(['PV-Budget Lower Level r=',num2str(R(ii))])

end
