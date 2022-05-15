function [partial_t,low,bp_low,pb_low,heat] = main_PV(w,Tau,Phi,R,dt,N,dx)

d_1 = d_1x(N,dx); d_2 = d_2x(N,dx); d_3 = d_3x(N,dx);

A2 = d_2x(N,dx); A2(1,:) = 1;

%% compare magnitude of adv. vs. heat. term in the dry lower-PV Equation:
%% d_t (phi_xx-tau_xx-tau)-phi_xxx+tau_xxx-phi_x = (1-r(w))*w.

phi = Phi(:,end); phi(1) = 0;
phi_p = Phi(:,end-1); phi_p(1)=0;
tau = Tau(:,end); tau(1) = 0;
tau_p = Tau(:,end-1); tau_p(1)=0;


phi = A2\phi; phi_p = A2\phi_p;

tau = A2\tau; tau_p = A2\tau_p;


% define r-factor based on w

r = r_factor(w,R);

% calculate adv. and heat. terms

q = d_2*phi-d_2*tau+tau;
q_p = d_2*phi_p-d_2*tau_p+tau_p;

partial_t = 1/dt*(q-q_p)-mean(r.*w);

low = d_3*tau-d_3*phi-d_1*phi;

bp_low = d_3*tau-d_3*phi;

pb_low = -d_1*phi;

heat = +(1-r).*w;



end