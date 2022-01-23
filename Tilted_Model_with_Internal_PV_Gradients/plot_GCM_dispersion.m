% Compare the 2-layer model numerical predictions with drag and beta on finite domain
% against the GCM results. 

close all; clear;

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


% load results of 2-layer models

list = {'_mid','_mid_beta_drag','_mid_beta','_mid_drag','_galerkin','_galerkin_beta_drag'};
path = '/disk7/mkohl/Mode_Kohl_OGorman_beta/2D_3D_comparison/';

for tt = 1:length(list)
    
load(strcat(path,'DRV',list{tt},'.mat'));
load(strcat(path,'MM',list{tt},'.mat'));
end

% Plot sigma and b for the paper

hfig = figure(1);
pos = get(hfig,'position');
set(hfig,'position',pos.*[.5 1 2 1])

subplot(1,2,1)

plot(surf_temp,sigma_GCM,'linewidth',1.5); hold on;
plot(surf_temp,rescale_sigma_mid_new.*sigma_DRV_mid_beta_drag,'linewidth',1.5); hold on;
plot(surf_temp,rescale_sigma_mid_new.*sigma_MM_mid_beta_drag,'linewidth',1.5); hold on;
xlabel('Global-mean surface air temperature (K)'); ylabel('\sigma (day^{-1})')
title('(a) Growth Rate')
legend('GCM','Tilted','Untilted','Location','NorthWest'); legend boxoff
set(gca,'FontSize',12)
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')


subplot(1,2,2)
plot(surf_temp,b_mid/1000,'linewidth',1.5); hold on;
plot(surf_temp(6:end),rescale_b_mid_new(6:end).*b_DRV_mid_beta_drag(6:end),'linewidth',1.5); hold on;
plot(surf_temp,rescale_b_mid_new.*b_MM_mid_beta_drag,'linewidth',1.5); hold on;
xlabel('Global-mean surface air temperature (K)'); ylabel('b (km)')
title('(b) Half-Ascent Area')
legend('GCM','Tilted','Untilted','Location','NorthEast'); legend boxoff
set(gca,'FontSize',12)
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')

%saveas(gcf,'/net/halo/disk28/disk7/mkohl/Mode_Kohl_OGorman/DRV_Paper/Figures/GCM_dispersion','epsc')

close all;

% Make Plots - without drag, and beta

figure(1)
plot(surf_temp,sigma_GCM,'linewidth',1.5); hold on;
plot(surf_temp,rescale_sigma_mid_new.*sigma_DRV_mid,'linewidth',1.5); hold on;
plot(surf_temp,rescale_sigma_mid_new.*sigma_MM_mid,'linewidth',1.5); hold on;
xlabel('T_g (K)')
ylabel('\sigma (day^{-1})')
title('Growthrate without beta and drag')
legend('GCM','DRV','MM','Location','NorthWest')
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);

figure(2)
plot(surf_temp,b_mid/1000,'linewidth',1.5); hold on;
plot(surf_temp,rescale_b_mid_new.*b_DRV_mid,'linewidth',1.5); hold on;
plot(surf_temp,rescale_b_mid_new.*b_MM_mid,'linewidth',1.5); hold on;
xlabel('T_g (K)')
ylabel('km')
title('b without beta and drag')
legend('GCM','DRV','MM','Location','NorthWest')
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);

% sigma and b with drag and beta

figure(3)
plot(surf_temp,sigma_GCM,'linewidth',1.5); hold on;
plot(surf_temp,rescale_sigma_mid_new.*sigma_DRV_mid_beta_drag,'linewidth',1.5); hold on;
plot(surf_temp,rescale_sigma_mid_new.*sigma_MM_mid_beta_drag,'linewidth',1.5); hold on;
xlabel('T_g (K)')
ylabel('\sigma (day^{-1})')
title('Growthrate with beta and drag')
legend('GCM','DRV','MM','Location','NorthWest')
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);

figure(4)
plot(surf_temp,b_mid/1000,'linewidth',1.5); hold on;
plot(surf_temp,rescale_b_mid_new.*b_DRV_mid_beta_drag,'linewidth',1.5); hold on;
plot(surf_temp,rescale_b_mid_new.*b_MM_mid_beta_drag,'linewidth',1.5); hold on;
xlabel('T_g (K)')
ylabel('km')
title('b with beta and drag')
legend('GCM','DRV','MM','Location','NorthWest')
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);

% just drag

figure(5)
plot(surf_temp,sigma_GCM,'linewidth',1.5); hold on;
plot(surf_temp,rescale_sigma_mid_new.*sigma_DRV_mid_drag,'linewidth',1.5); hold on;
plot(surf_temp,rescale_sigma_mid_new.*sigma_MM_mid_drag,'linewidth',1.5); hold on;
xlabel('T_g (K)')
ylabel('\sigma (day^{-1})')
title('Growthrate with drag')
legend('GCM','DRV','MM','Location','NorthWest')
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);

figure(6)
plot(surf_temp,b_mid/1000,'linewidth',1.5); hold on;
plot(surf_temp,rescale_b_mid_new.*b_DRV_mid_drag,'linewidth',1.5); hold on;
plot(surf_temp,rescale_b_mid_new.*b_MM_mid_drag,'linewidth',1.5); hold on;
xlabel('T_g (K)')
ylabel('km')
title('b with drag')
legend('GCM','DRV','MM','Location','NorthWest')
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);

% just beta

figure(7)
plot(surf_temp,sigma_GCM,'linewidth',1.5); hold on;
plot(surf_temp,rescale_sigma_mid_new.*sigma_DRV_mid_beta,'linewidth',1.5); hold on;
plot(surf_temp,rescale_sigma_mid_new.*sigma_MM_mid_beta,'linewidth',1.5); hold on;
xlabel('T_g (K)')
ylabel('\sigma (day^{-1})')
title('Growthrate with beta')
legend('GCM','DRV','MM','Location','NorthWest')
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);

figure(8)
plot(surf_temp,b_mid/1000,'linewidth',1.5); hold on;
plot(surf_temp,rescale_b_mid_new.*b_DRV_mid_beta,'linewidth',1.5); hold on;
plot(surf_temp,rescale_b_mid_new.*b_MM_mid_beta,'linewidth',1.5); hold on;
xlabel('T_g (K)')
ylabel('km')
title('b with beta')
legend('GCM','DRV','MM','Location','NorthWest')
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);

% Galerkin without beta and drag

figure(9)
plot(surf_temp,sigma_GCM,'linewidth',1.5); hold on;
plot(surf_temp,rescale_sigma_galerkin.*sigma_DRV_galerkin,'linewidth',1.5); hold on;
plot(surf_temp,rescale_sigma_galerkin.*sigma_MM_galerkin,'linewidth',1.5); hold on;
xlabel('T_g (K)')
ylabel('\sigma (day^{-1})')
title('Growthrate without beta and drag: Galerkin')
legend('GCM','DRV','MM','Location','NorthWest')
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);

figure(10)
plot(surf_temp,b_mid/1000,'linewidth',1.5); hold on;
plot(surf_temp,rescale_b_galerkin.*b_DRV_galerkin,'linewidth',1.5); hold on;
plot(surf_temp,rescale_b_galerkin.*b_MM_galerkin,'linewidth',1.5); hold on;
xlabel('T_g (K)')
ylabel('km')
title('b without beta and drag: Galerkin')
legend('GCM','DRV','MM','Location','NorthWest')
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);

figure(11)
plot(surf_temp,sigma_GCM,'linewidth',1.5); hold on;
plot(surf_temp,rescale_sigma_galerkin.*sigma_DRV_galerkin_beta_drag,'linewidth',1.5); hold on;
plot(surf_temp,rescale_sigma_galerkin.*sigma_MM_galerkin_beta_drag,'linewidth',1.5); hold on;
xlabel('T_g (K)')
ylabel('\sigma (day^{-1})')
title('Growthrate: Galerkin')
legend('GCM','DRV','MM','Location','NorthWest')
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);

figure(12)
plot(surf_temp,b_mid/1000,'linewidth',1.5); hold on;
plot(surf_temp,rescale_b_galerkin.*b_DRV_galerkin,'linewidth',1.5); hold on;
plot(surf_temp,rescale_b_galerkin.*b_MM_galerkin,'linewidth',1.5); hold on;
xlabel('T_g (K)')
ylabel('km')
title('b: Galerkin')
legend('GCM','DRV','MM','Location','NorthWest')
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);

