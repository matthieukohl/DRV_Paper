function [partial_t,low,bp_low,pb_low,heat] = main_PV_Budget(w,tau,tau_p,phi,phi_p,R,dt,N,dx,beta)

d_1 = d_1x(N,dx); d_2 = d_2x(N,dx); d_3 = d_3x(N,dx);

% Calculate Terms in the lower PV-Budget: PV_2 = (phi-tau)_xx + tau + h
% Eulerian -------------- + (u*q_x)                 + (v*q_y) 
% d_t (phi_xx-tau_xx+tau) + (tau_xxx-phi_xxx-tau_x) + (tau+h)_y*(phi-tau)_x
% = Heating 
% = (1-r(w))*w

% Note that v'*q_y = (tau+h)_y *(phi-tau)_x = 0 by construction

% define r-factor based on w

r = r_factor(w,R);

% calculate adv. and heat. terms

q = d_2*phi-d_2*tau+tau;
q_p = d_2*phi_p-d_2*tau_p+tau_p;

partial_t = 1/dt*(q-q_p)-mean(r.*w);

bp_low = d_3*tau-d_3*phi-d_1*tau;

pb_low = beta*d_1*(phi-tau);

low = bp_low + pb_low;

heat = +(1-r).*w;

end