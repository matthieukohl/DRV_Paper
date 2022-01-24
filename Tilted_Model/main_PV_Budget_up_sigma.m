function [partial_t,up,bp_up,pb_up,heat] = main_PV_Budget_up_sigma(w,tau,phi,R,sigma,N,dx)

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

q = d_2*phi+d_2*tau-tau;

partial_t = sigma*q+mean(r.*w);

bp_up = d_3*tau+d_3*phi-d_1*tau;

pb_up = zeros(N,1);

up = bp_up + pb_up;

heat = -(1-r).*w;

end