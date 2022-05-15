% Mode or DRV plot

close all; clear;

r = {'r001','r01'};

for ii = 1:2

load(['mode_or_drv/mode_or_drv_h1_vs_h2_fixed_',r{ii},'_N100_L8pi_t200.mat'])

% set NaN for all negative growth rates

mode_or_drv(sigma<0) = NaN;

% define the pv gradients based on h1 and h2

q1y = 1-h1; q2y = -1+h2;

[Q1y,Q2y] = meshgrid(q1y,q2y);


x0=10; y0=10; width=1000; height=400; 
set(gcf,'position',[x0,y0,width,height])

% Plot mode vs DRV
subplot(1,2,ii)
contourf(Q1y,Q2y,mode_or_drv,[-1 1],'Edgecolor','none'); hold on
contour(Q1y,Q2y,sigma,'k--','ShowText','on')
colormap redblue
colorbar
%legend('Location','eastoutside')
set(gca,'Ydir','reverse')
xlabel('q_{1y}'); ylabel('q_{2y}')
title(['R=',num2str(R)])
set(gca,'FontSize',12)
caxis([-5 5])
%brighten(0.8)
end


% try to make the plot as scatter 
x_drv = [];
y_drv = [];
x_mode = [];
y_mode = [];

for ii  = 1:length(h1)
    for jj = 1:length(h2)
        
       if mode_or_drv(ii,jj)==1
           x_drv = cat(1,x_drv, q1y(ii));
           y_drv = cat(1,y_drv,q2y(jj));
       elseif mode_or_drv(ii,jj) == -1
           x_mode = cat(1,x_mode,q1y(ii));
           y_mode = cat(1,y_mode,q2y(jj));
       end
    end
end

figure(2)
scatter(x_mode,y_mode,'b'); hold on;
scatter(x_drv,y_drv,'r');
%contour(Q1y,Q2y,sigma,'k--','ShowText','on')
set(gca,'Ydir','reverse')
