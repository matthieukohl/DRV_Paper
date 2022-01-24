function [rescale_sigma_galerkin,rescale_b_galerkin,rescale_sigma_av,rescale_b_av,...
 rescale_sigma_omega,rescale_b_omega,rescale_sigma_500,rescale_b_500,...
 rescale_sigma_400,rescale_b_400,rescale_sigma_mid,rescale_b_mid,rescale_sigma_mid_av,rescale_b_mid_av,...
 rescale_sigma_wmax,rescale_b_wmax,rescale_sigma_LH,rescale_b_LH,rescale_sigma_mid_new,rescale_b_mid_new] = ...
 rescale_2D_3D(u_galerkin,u_mean,u_omega,u_500,u_400,u_mid,u_wmax,u_LH,u_mid_new,...
 f_eff,delta_p,delta_p_500,delta_p_400,delta_p_mid,delta_p_wmax,delta_p_LH,delta_p_mid_new,...
S_galerkin,S_mean,S_omega,S_500,S_400,S_mid,S_mid_av,S_wmax,S_LH,S_mid_new)

% Calculate the rescale factors 

day = 24*3600; km = 1000;

% Galerkin weighting

%rescale_sigma_galerkin = u_galerkin.*f_eff./(delta_p.*sqrt(-S_galerkin))*day;

%rescale_b_galerkin = sqrt(-S_galerkin).*delta_p./f_eff*1/km;


rescale_sigma_galerkin = u_galerkin.*f_eff./((delta_p/2).*sqrt(-S_galerkin))*day;

rescale_b_galerkin = sqrt(-S_galerkin).*(delta_p/2)./f_eff*1/km;

% Mean

rescale_sigma_av = u_mean.*f_eff./(delta_p.*sqrt(-S_mean))*day;

rescale_b_av = sqrt(-S_mean).*delta_p./f_eff*1/km;

% 500hPA

rescale_sigma_500 = u_500.*f_eff./(delta_p_500.*sqrt(S_500/2))*day;

rescale_b_500 = sqrt(S_500/2).*delta_p_500./f_eff*1/km;

% omega-weighting

rescale_sigma_omega = u_omega.*f_eff./(delta_p.*sqrt(S_omega))*day;

rescale_b_omega = sqrt(S_omega).*delta_p./f_eff*1/km;

% midlevel - weighting

rescale_sigma_mid = u_mid.*f_eff./(delta_p_mid.*sqrt(S_mid/2))*day;

rescale_b_mid = sqrt(S_mid/2).*delta_p_mid./f_eff*1/km;

% midlevel-average

rescale_sigma_mid_av = u_mid.*f_eff./(delta_p_mid.*sqrt(S_mid_av/2))*day;

rescale_b_mid_av = sqrt(S_mid_av/2).*delta_p_mid./f_eff*1/km;

% midlevel new - weighting

rescale_sigma_mid_new = u_mid_new.*f_eff./(delta_p_mid_new.*sqrt(S_mid_new/2))*day;

rescale_b_mid_new = sqrt(S_mid_new/2).*delta_p_mid_new./f_eff*1/km;

% Latent heating weighting

rescale_sigma_LH = u_LH.*f_eff./(delta_p_LH.*sqrt(S_LH/2))*day;

rescale_b_LH = sqrt(S_LH/2).*delta_p_LH./f_eff*1/km;

% Galerkin with sine weighting based on wmax

rescale_sigma_wmax = u_wmax.*f_eff./((delta_p_wmax/2).*sqrt(-S_wmax))*day;

rescale_b_wmax = sqrt(-S_wmax).*(delta_p_wmax/2)./f_eff*1/km;

% H topped at 400hPa 

rescale_sigma_400 = u_400.*f_eff./(delta_p_400.*sqrt(S_400/2))*day;

rescale_b_400 = sqrt(S_400/2).*delta_p_400./f_eff*1/km;

end

