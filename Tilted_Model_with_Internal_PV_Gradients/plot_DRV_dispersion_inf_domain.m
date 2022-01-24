% Compare the root-solving solutions for sigma and b of the DRV dispersion 
% equations on an infinite domain against the GCM results.

clear; close all;

% load the r values inferred from the GCM

load('2D_3D_data.mat','R_mid_new');

r = flip(R_mid_new); %clear('R_mid_new')

r(1) = 1e-3; 
r = r(1:8);
r = r';


% Calculate sigma and b from dispersion relations on an infinite domain

% Define Solver Options

options = optimset('TolX',1e-12,'TolFun',1e-12);
 
% Define Initial Guess

%x0(1) = 1.5; % first solution
%x0(1) = 1.6; % second solution

%x0(2) = 0.2; % first solution
%x0(2) = 0.1; % second solution



% Pre-define values

[s,b,fval] = deal(zeros(1,length(r)));

% asymptotic predictions: sigma = c0, b = pi/2*sqrt(r)+(1+c0)*r

c0 = 1/2*(1+sqrt(5)); c1 = pi*(1+c0)*(c0^2-1)/(2*(1-2*c0));

d2 = c1+pi*(1+2*c0^2)/4;

s_asym = c0 + c1*sqrt(r)  ;
b_asym = pi/2*sqrt(r)+(1+c0)*r + d2*r.^3/2;


% s_asym = 1.545 -1.852*sqrt(r) + 1.037*r;
% 
% b_asym = 28.42*sqrt(r) -161.2*r + 236*r.^(3/2);



 
 for i = 1:length(r)

% linear interpolation as initial guess

x0(1) = s_asym(i); 
x0(2) = b_asym(i);

if i>6
x0(1) = 2*s(i-1)-s(i-2); 
x0(2) = 2*b(i-1)-b(i-2); 
end

s0(i) = x0(1); b0(i) = x0(2);


% Solve the tan-equations for growth rate s and half-ascent area b

fun = @(x)residual_fn(x,r(i));
[x,fval(i)] = fminsearch(fun,x0,options);

s(i) = x(1);
b(i) = x(2);


end
 
% Calculate and plot the equation residuals
[k1,k2] = Wavenumbers(s,r);

[Eq1,Eq2] = tan_equations(k1,k2,s,b,r);

figure(1)
set(gcf,'position',[20 50 1250 600])
subplot(1,2,1)
loglog(r,abs(Eq1))
xlabel('r');ylabel('residual')
set(gca,'FontSize',14)
title('Tan-Equation1')

subplot(1,2,2)

loglog(r,abs(Eq2))
xlabel('r');ylabel('residual')
set(gca,'FontSize',14)
title('Tan-Equation2')

% flip sigma and b to go from cold to warm

s = flip(s); b = flip(b);

% Load in the rescale factors to convert to dimensional units
surf_temp = [270.0311,277.5996,280.5327,283.1391,285.4284,287.6549,290.8506...
293.6622, 296.1045, 298.1490, 300.0797, 303.7376, 306.5808, ...
310.7297, 316.1383];

load('/net/halo/disk28/disk7/mkohl/Mode_Kohl_OGorman_beta/data.mat')

[rescale_sigma_galerkin,rescale_b_galerkin,rescale_sigma_av,rescale_b_av,...
 rescale_sigma_omega,rescale_b_omega,rescale_sigma_500,rescale_b_500,...
 rescale_sigma_400,rescale_b_400,rescale_sigma_mid,rescale_b_mid,rescale_sigma_mid_av,rescale_b_mid_av,...
 rescale_sigma_wmax,rescale_b_wmax,rescale_sigma_LH,rescale_b_LH,rescale_sigma_mid_new,rescale_b_mid_new] = ...
 rescale_2D_3D(u_galerkin,u_mean,u_omega,u_500,u_400,u_mid,u_wmax,u_LH,u_mid_new,...
 f_eff,delta_p,delta_p_500,delta_p_400,delta_p_mid,delta_p_wmax,delta_p_LH,delta_p_mid_new,...
S_galerkin,S_mean,S_omega,S_500,S_400,S_mid,S_mid_av,S_wmax,S_LH,S_mid_new);


hfig = figure(1);
%pos = get(hfig,'position');
%set(hfig,'position',pos.*[.5 1 2 1])

subplot(1,2,1)

plot(surf_temp,sigma_GCM,'linewidth',1.5); hold on;
plot(surf_temp(8:end),rescale_sigma_mid_new(8:end).*s','linewidth',1.5); hold on;
xlabel('Global-mean surface air temperature (K)'); ylabel('\sigma (day^{-1})')
title('(a) Growth Rate')
legend('GCM','Tilted','Location','NorthWest'); legend boxoff
set(gca,'FontSize',12)
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')


subplot(1,2,2)
plot(surf_temp,b_mid/1000,'linewidth',1.5); hold on;
plot(surf_temp(8:end),rescale_b_mid_new(8:end).*b','linewidth',1.5); hold on;
xlabel('Global-mean surface air temperature (K)'); ylabel('b (km)')
title('(b) Half-Ascent Area')
legend('GCM','Tilted','Untilted','Location','NorthEast'); legend boxoff
set(gca,'FontSize',12)
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')


