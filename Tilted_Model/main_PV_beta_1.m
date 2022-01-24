function [partial_t,low,bp_low,pb_low,heat,diffusion] = main_PV_beta_1(w,tau,tau_p,phi,phi_p,R,beta,nu,dt,N,dx)

d_1 = d_1x(N,dx); d_2 = d_2x(N,dx); d_3 = d_3x(N,dx); d_4 = d_4x(N,dx);

A2 = d_2x(N,dx); A2(1,:) = 1;

% d_1 = d_1x_4(N,dx); d_2 = d_2x_4(N,dx); d_3 = d_3x_4(N,dx); d_4 = d_4x_4(N,dx);
% 
% A2 = d_2x_4(N,dx); A2(1,:) = 1;

%% compare magnitude of adv. vs. heat. term in the dry lower-PV Equation:
%% d_t (phi_xx-tau_xx-tau)-phi_xxx+tau_xxx-phi_x = (1-r(w))*w.

% define r-factor based on w

r = r_factor(w,R);

% calculate adv. and heat. terms

q = d_2*phi-d_2*tau+tau;
q_p = d_2*phi_p-d_2*tau_p+tau_p;

partial_t = 1/dt*(q-q_p)-mean(r.*w);

low = d_3*tau-d_3*phi-d_1*phi+beta*(d_1*phi-d_1*tau);

bp_low = d_3*tau-d_3*phi-d_1*tau;

pb_low = -(d_1*phi-d_1*tau)+beta*(d_1*phi-d_1*tau);

heat = +(1-r).*w;

diffusion = nu*d_4*(phi-tau);

end