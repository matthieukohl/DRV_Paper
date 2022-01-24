% Mode or DRV plot

close all; clear;

%load('mode_or_drv/mode_or_drv_L16pi_t400.mat')

%load('mode_or_drv/mode_or_drv_N400_L32pi_t400.mat')

%load('mode_or_drv/mode_or_drv_N200_L8pi_t200_new.mat')

load('mode_or_drv/mode_or_drv_N400_L16pi_t200_new.mat')

% mode_or_drv(1:end-1,3) = NaN;
% mode_or_drv(1:end-1,4) = NaN;
%mode_or_drv(end,11:end) = NaN;
%mode_or_drv(end,1:10) = -1;
% 
% index = find(h==1);
% sigma(end,index:end) = 0;
% mode_or_drv(end,index:end) = NaN;
% mode_or_drv(end,index-2) = 1;
% mode_or_drv(1,6) = 1;


% sigma(end,19:end) = 0;
% mode_or_drv(sigma<=0) = NaN;
% b(sigma<=0) = NaN;
% sigma(sigma<=0) = NaN;

qy = 1-h;

[RR,hh] = meshgrid(qy,R);

figure(1)
contourf(hh,RR,mode_or_drv,[1 -1]); hold on
contour(hh,RR,sigma,'k--','ShowText','on')
set(gca,'Ydir','reverse')
xlabel('r'); ylabel('qy_{up}')
title('Phase Diagram')
set(gca,'FontSize',12)
%saveas(gcf,'/net/halo/disk28/disk7/mkohl/Mode_Kohl_OGorman/DRV_Paper/Figures/mode_or_drv','epsc')


drv = [-1,1];

figure(2)
contourf(hh,RR,sigma); hold on;
%contour(hh,RR,mode_or_drv,drv,'k--','ShowText','on')
colorbar
caxis([-max(sigma(:)) max(sigma(:))])
colormap(redblue)
set(gca,'Ydir','reverse')
xlabel('r'); ylabel('qy_{up}')
title('sigma')

figure(3)
contourf(hh,RR,b); colorbar
xlabel('r'); ylabel('qy_{up}')
set(gca,'Ydir','reverse')
title('b')

% for ii = 1:length(R)
%    figure;
%    plot(squeeze(ww(ii,4,:)))
% end