% Solve the Mode-Equation as an initial value problem:
% d_tt ((rw)_xx-4*w) = (rw)_xxxx+4*w_xx
% d_tt B*w = A*w -> d_tt w = B^(-1)*A*w

clear;
close all;

% Define Grid

N = 200; 
L = 8*pi; 
x = linspace(0,L,N+1);
dx = L/(N);

% Define time-step

t_final = 40;
tN = 1600;
dt = t_final/tN;
time = linspace(0,t_final,tN);

% define Wave-numbers for further analysis

k0 =  2*pi/L;
k = [0:N/2,-N/2+1:-1]*k0;

% define r-factor

R = [1.0;0.4;0.10;0.04;0.01;-0.1]; 

% define variables

sigma = zeros(length(R),1); lambda = zeros(length(R),1); 

w_final = zeros(N+1,length(R));

[partial_t,adv_up,adv_low,adv_bp_up,adv_bp_low,adv_pb_up,adv_pb_low,heating,...
budget_w,budget_w_mode,budget_PV] = deal(zeros(length(R),1));

[w_equation,w_equation_mode,PV_equation,...
    phi_equation,tau_equation] = deal(zeros(N,length(R))); 

% Initialize Solution

w_new = sin(x)+sin(2*x)+sin(4*x)+sin(6*x)+sin(8*x)+sin(10*x); 
w_new = w_new(1:end-1);


for ii = 6:6
ii

% Initialize Values

[w_test,r_test,phi_test,tau_test] = deal(zeros(N,3));

% Initial values    
w = sin(x)+sin(2*x)+sin(4*x)+sin(6*x)+sin(8*x)+sin(10*x); w = w(1:end-1);
w_t = zeros(1,N);  
w = [w,w_t]';  

%load('vortex.mat','w');

%Initial Values for Phi

phi = zeros(1,N);%sin(x)+sin(3*x)+sin(5*x)+sin(7*x)+sin(9*x)+sin(10*x); phi = phi(1:end-1);
tau = zeros(1,N);%cos(x)+sin(3*x)+cos(5*x)+cos(7*x)+sin(9*x)+cos(10*x); tau = tau(1:end-1);
Phi = [phi,tau]';

%Propagator for tau, phi fields

A1 = [zeros(N,N),-d_1x(N,dx);-d_1x(N,dx),zeros(N,N)];
A2 = d_2x(N,dx); A2(1,:) = 1;
%A3 = [zeros(N,N),-d_3x(N,dx);-d_3x(N,dx),zeros(N,N)];
%A_2 = d_2x(N,dx); A_2(1,:) = 1;
%A2 = [A_2,zeros(N,N);zeros(N,N),A_2]; A2_inv = inv(A2);

Propagator = diag(ones(1,2*N))-dt*A1; 
%Propagator = diag(ones(1,2*N))-dt*A2_inv*A3; 

count = 0;
index = 0;

%[t,y] = ode15s(@(t,y) Evolution(t,y,R(ii),N,dx),[0 t_final],w);
    
for tt = 1:tN
tt
count = count +1;   
% Define r-factor

r = r_factor(w(1:N),R(ii));

% Define Propagator C = [0,1;0,B^(-1)*A]

A = A_matrix_correct(r,N,dx); B = B_matrix(r,N,dx); 

%A = 100*d_2x(N,dx); 

% C = [zeros(N,N),diag(ones(1,N));inv(B)*A,zeros(N,N)];
% 
% D = diag(ones(1,2*N))-dt*C; 

% alternative time differencing

C = [zeros(N,N),inv(B);A,zeros(N,N)];

D = diag(ones(1,2*N))-dt*C; 

% Propagate solution forward using implicit euler and rescale

% w = inv(D)*w;
% 
% w = w/norm(w);

% Propage solution forward using SVD. 

% Propag = expm(C*1); 
% [U,S,V] = svd(Propag);
% w = U(1:N,1);
% w = w/norm(w);

% Use weighting 

%w = 0.5*inv(D)*w+0.5*w;

w = inv(D)*w;

plot(w(1:N)); pause(0.01);

%w = w/norm(w);

% Use ODE-45
% 
% [t,w] = ode45(@Propagtor(,[time time+dt],w);
% time = time+dt;

% Propagate tau and phi fields

c = [zeros(N,1);-2*w(1:N)];
% %c = A2_inv*[zeros(N,1);-2*w(1:N)];

Phi = inv(Propagator)*(Phi+dt*c);

%Plot evolution of w
%if mod(tt,200)==0
% plot(x(1:N),w(1:N)); 
% pause(0.1)
% hold off;
%end

%

Stream1 = Phi(1:N); Stream1(1) = 0;
Stream2 = Phi(N+1:end); Stream2(1) = 0;
phi(tt,:) = A2\Stream1;
tau(tt,:) = A2\Stream2;
W(tt,:) = w(1:N);

% Store last three w-values to check validity of equation

if tt == tN-2 || tt == tN-1 || tt == tN
   
   index = index+1;
   w_test(:,index) = w(1:N);
   r_test(:,index) = r_factor(w(1:N),R(ii));
   Phi_test(:,index) = Phi;
   
   
   Stream1 = Phi(1:N); Stream1(1) = 0;
   Stream2 = Phi(N+1:end); Stream2(1) = 0;
   phi_test(:,index) = A2\Stream1;
   tau_test(:,index) = A2\Stream2;
  
%    if index == 2 || index == 3
%    growthrate = 1/dt*log(mean(w_test(:,index)./w_test(:,index-1)));
%    
%    Stream1 = Phi(1:N); Stream1(1) = 0;
%    Stream2 = Phi(N+1:end); Stream2(1) = -1/(2*growthrate)*mean(r_test(:,index).*w_test(:,index));
%    
%    phi_test(:,index) = A2\Stream1;
%    tau_test(:,index) = A2\Stream2;
%    
%    end
   
   
end

end

%w_new = w(1:N)';

% Store final w-fields and the asymmetry factor

w_final(:,ii) = [w(1:N);w(1)];
lambda(ii) = asymmetry(w_final(:,ii)/norm(w_final(:,ii)));

% Calculate growthrate assuming exponential increase at late time

%sigma(ii) = 1/dt*log(mean(w_test(:,3)./w_test(:,2)));

sigma(ii) = growthrate(phi,time);

% Calculate residual of d_tt ((rw)_xx-w) = ((rw)_4x+w_2x)
% using backward Euler

RHS  = A_matrix_correct(r_test(:,3),N,dx)*w_test(:,3);
LHS = B_matrix(r_test(:,3),N,dx)*w_test(:,3);
LHSm1 = B_matrix(r_test(:,2),N,dx)*w_test(:,2);
LHSm2 = B_matrix(r_test(:,1),N,dx)*w_test(:,1);

% How well does the w-equation close?

w_equation(:,ii) = 1/(norm(w_test(:,3)))*((LHS+LHSm2-2*LHSm1)/dt^2-RHS);

w_equation_mode(:,ii) = 1/norm(w_test(:,3))*(sigma(ii)^2*B_matrix(r_test(:,3),N,dx)*w_test(:,3)...
    -A_matrix_correct(r_test(:,3),N,dx)*w_test(:,3));

tau_equation(:,ii) = 1/norm(tau_test(:,3))*(1/dt*(A2*(tau_test(:,3)-tau_test(:,2)))+d_3x(N,dx)*phi_test(:,3)+2*w_test(:,3));

phi_equation(:,ii) = 1/norm(phi_test(:,3))*(1/dt*(A2*(phi_test(:,3)-phi_test(:,2)))+d_3x(N,dx)*tau_test(:,3));

Phi_equation(:,ii) = 1/norm(Phi_test(:,3))*(1/dt*(Phi_test(:,3)-Phi_test(:,2))-A1*Phi_test(:,3)+[zeros(N,1);2*w_test(:,3)]);

budget_w(ii) = rms(w_equation(:,ii)); budget_w_mode(ii) = rms(w_equation_mode(:,ii));

% Calculate PV_Budget

%PV2 = PV_spec(w(1:N)'/norm(w(1:N)'),sigma(ii),L,N);

%q = PV(w(1:N)/norm(w(1:N)),sigma(ii),N,dx);

%[partial_t_term,up,low,bp_up,bp_low,pb_up,pb_low,heat] = Diagnostic_shape(w(1:N)/norm(w(1:N)),R(ii),sigma(ii),L,N);

%[partial_t_term,low,bp_low,pb_low,heat] = q_Budget(w(1:N),R(ii),sigma(ii),N,dx);

[partial_t_term,low,bp_low,pb_low,heat] = PV_Budget(w_test(:,3),tau_test(:,3),tau_test(:,2),phi_test(:,3),phi_test(:,2),R(ii),dt,N,dx);

adv_low(ii) = rms(low); adv_bp_low(ii) = rms(bp_low); adv_pb_low(ii) = rms(pb_low);

partial_t(ii) = rms(partial_t_term); heating(ii) = rms(heat);

% How well does the lower layer PV-equation close?

PV_equation(:,ii) = (partial_t_term+pb_low+bp_low-heat);

budget_PV(ii) = rms(PV_equation(:,ii));


% Plot Various fields

% Plot adv. - heat. profiles

% close all;
% 
figure(24)
plot(x,-1/norm(partial_t_term)*[heat;heat(1)]); hold on;
plot(x,1/norm(partial_t_term)*[partial_t_term;partial_t_term(1)]); hold on
plot(x,1/norm(partial_t_term)*[bp_low;bp_low(1)]); hold on;
plot(x,1/norm(partial_t_term)*[pb_low;pb_low(1)]); hold on
plot(x,1/norm(partial_t_term)*[PV_equation(:,ii);PV_equation(1,ii)])
plot(x,1/norm(partial_t_term)*w_final(:,ii)); hold on;
xlabel('x'); ylabel('Magnitude')
legend('heat','Partial t','uqx','vqy','Residual','w');
title(['PV-Budget Lower Level r=',num2str(R(ii))])
%xlim([0 4*pi])
% savefig(['/net/aimsir/archive1/mkohl/Mode_Zurita_Gotor/figures/initial_value/PV_Budget_rand_r',num2str(R(ii)),'.fig'])
% close all;
% 
% figure(26)
% plot(x,w_final(:,ii))
% xlabel('x'); ylabel('w')
% title(['w-profile r=',num2str(R(ii))])
% savefig(['/net/aimsir/archive1/mkohl/Mode_Zurita_Gotor/figures/initial_value/w_rand_r',num2str(R(ii)),'.fig'])
% close all;



% 
% % % 
% figure(27)
% 
% plot(x,[q;q(1)]); hold on;
% xlabel('x'); ylabel('advection')
% legend('PV');
% title(['PV Lower Level r=',num2str(R(ii))])
% savefig(['/net/aimsir/archive1/mkohl/Mode_Zurita_Gotor/figures/initial_value/alternative/fields/PV_rand_r',num2str(R(ii)),'.fig'])
% close all;
% % % 
% % figure(26)
% % plot(x,w_final)
% % xlabel('x'); ylabel('w')
% % title(['w-profile r=',num2str(R(ii))])
% % savefig(['/net/aimsir/archive1/mkohl/Mode_Zurita_Gotor/figures/initial_value/alternative/fields/w_rand_r',num2str(R(tt)),'.fig'])
% % close all;
% % % 
% figure(28)
% plot(x,[w_equation(:,ii);w_equation(1,ii)])
% xlabel('x'); ylabel('Sum w-equation')
% title(['w-equation r=',num2str(R(ii))])
% savefig(['/net/aimsir/archive1/mkohl/Mode_Zurita_Gotor/figures/initial_value/alternative/fields/w_equation_rand_r',num2str(R(ii)),'.fig'])
% close all;
% % % 
% figure(29)
% plot(x,[PV_equation(:,ii);PV_equation(1,ii)])
% xlabel('x'); ylabel('Sum PV-equation')
% title(['PV-equation r=',num2str(R(ii))])
% savefig(['/net/aimsir/archive1/mkohl/Mode_Zurita_Gotor/figures/initial_value/alternative/fields/PV_equation_16pi_r',num2str(R(ii)),'.fig'])
% close all;


end

% f = w_final(1:end-1,1);
% f_t = fft(f);
% Power = abs(f_t).^2;
% a = find(max(Power)==Power);
% k(a(1))


% reduction_factor = [0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.6,0.8,1.0];
% 
% 
% lambda_omega_mode = [0.9283,0.9173,0.8996,0.8313,0.7820,0.7343,0.6871,...
%  0.5674,0.5187,0.4994];
% 
% lambda_omegaQG_mode = [0.9193,0.9053,0.8903,0.8200,0.7674,0.7277,0.6756,...
%  0.5706,0.5192,0.5043];
% % 
% % lambda_theory_mode = [0.5002, 0.5094, 0.5271, 0.5452, 0.5709, 0.6003, 0.6372, 0.6810,...
% % 0.7448, 0.7547, 0.7733, 0.7882, 0.8063, 0.8169, 0.8296, 0.8520, 0.8755, 0.9007,0.9317, 0.9530];
% % 
% % % 
% figure(1)
% semilogx(reduction_factor,lambda_omega_mode,'Color','red'); hold on;
% semilogx(R,lambda,'Color','blue'); hold on;
% xlabel('r')
% ylabel('Lambda')
% title('Asymmetry Modal Equation')
% legend('GCM','1D Equation')
% savefig('/net/aimsir/archive1/mkohl/Mode_Zurita_Gotor/figures/initial_value/lambda_rand_r.fig')
% % 
% close all;
% % 
% figure(2)
% plot(R,sigma/2)
% xlabel('r')
% ylabel('Lambda')
% title('Growthrate')
% savefig('/net/aimsir/archive1/mkohl/Mode_Zurita_Gotor/figures/initial_value/Growthrate_rand.fig')
% close all;
% % % 
% figure(3)
% semilogx(R,partial_t); hold on;
% semilogx(R,adv_bp_low); hold on;
% semilogx(R,adv_pb_low); hold on;
% semilogx(R,heating); 
% xlabel('r'); ylabel('rms')
% title('Magnitude of advection terms lower layer vs. heating term')
% legend('Eulerian','uqx_{low}','vqy_{low}','heat');
% savefig('/net/aimsir/archive1/mkohl/Mode_Zurita_Gotor/figures/initial_value/Norm_rand.fig')

% 
% 
% % % % 
% % % % 
% % % % % figure(2)
% % % % % for ii = 1:length(R)
% % % % %     
% % % % % plot(x,w_final(:,ii)); 
% % % % % xlabel('x'); ylabel('w')
% % % % % title(['w-field r=',num2str(R(ii))])
% % % % % pause(2.0)
% % % % % 
% % % % % end
% % % % 
% 
% 
% close all;
% 
% v = VideoWriter('/net/aimsir/archive1/mkohl/Mode_Zurita_Gotor/videos/initial_value/alternative/w_16pi_N200_16000.avi');
% v.FrameRate = 1;
%  
% open(v);
% 
% for tt = 1:length(R)
%     
% plot(x,w_final(:,tt)); 
% xlabel('x')
% ylabel('w')
% title(['r=',num2str(R(tt))]);
% xlim([0 L])
% 
% 
% frame = getframe(gcf);
% myVideo.FrameRate = 1.0;
% writeVideo(v,frame);
%   
% end
%     
% 
% close(v);
% % 
% % 
% % % 
