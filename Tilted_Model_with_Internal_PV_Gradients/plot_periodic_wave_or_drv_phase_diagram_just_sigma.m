% Mode or DRV plot

close all; clear;

x0=10; y0=10; width= 8/7*1100; height=8/7*300; 
set(gcf,'position',[x0,y0,width,height])


load('mode_or_drv/mode_or_drv_N200_L8pi_t200_new.mat')
%load('mode_or_drv/mode_or_drv_N50_L8pi_t200.mat')

% set all the stable modes to NaN values

mode_or_drv(sigma<0.09) = 2;
%mode_or_drv(sigma<0.09) = 2;


% if zero-local max set to DRV (lies at the boundary)

for ii = 1:length(R)
    for jj = 1:length(h)
ind = islocalmax(ww(ii,jj,:)); 
ind(ww(ii,jj,:)<0) = 0; 

if sum(ind)==0
    mode_or_drv(ii,jj) = 1;
end
    end
end
      
qy = 1-h;

%[RR,hh] = meshgrid(-qy,R);

[RR,hh] = meshgrid(qy,R);

levels = [1.4,1.2,0.8,0.6,0.4,0.09];

%b(sigma<0.09) = NaN;

%contourfcmap(hh,RR,mode_or_drv,[1 2 -1]','EdgeColor','none'); hold on

subplot(1,3,1)

contourf(hh,RR,mode_or_drv,[-1 0 2],'EdgeColor','none'); hold on
contour(hh,RR,sigma,levels,'k--','ShowText','on'); hold on; 
%contour(hh,RR,b,'g:','ShowText','on');
%colormap redblue
%map = [184/255 215/255 240/255;1 152/255 145/255;  1  1 1]; colormap (map);
map = [200/255 197/255 252/255;255/255 198/255 199/255;  1  1 1]; colormap (map);
caxis([-1 2]); % this is needed so that three colors only span variable mode_or_drv
% index = 3;
% rgba = double(cat(2, h.FacePrims.ColorData))./255;
% col = rgba(1:3,:)';
% col(index,:)=ones(size(col(index,:))); % [1 1 1] is the RGB color code for white
% colormap(col)

cmp = colormap;
%colorbar

%contour(hh,RR,sigma,'k--','ShowText','on'); hold on; 
title('\rm q_{2y} = -q_{1y}')
%set(gca,'Ydir','reverse')
xlabel('r'); %ylabel('q_{1y}=-q_{2y}')
ylabel('q_{1y}')
set(gca,'FontSize',12)
%set(gca,'clim',[-1 1])
%caxis([-5 5])

annotation('textbox', [0.07, 0.98, 0, 0], 'string', '(a)','fontsize',12)



annotation('textbox', [0.24, 0.89, 0.2, 0], 'string', 'Periodic Wave','Edgecolor','none','fontsize',12)
annotation('textbox', [0.24, 0.73, 0, 0], 'string', 'DRV','fontsize',12)
annotation('textbox', [0.27, 0.42, 0, 0], 'string', 'Stable','fontsize',12)




% % plot root terms for test
% 
% n = 7;
% 
% term1 = 1-2*h-sigma(n,:).^2;
% term2 = 1-6*sigma(n,:).^2+4*h.*sigma(n,:).^2+sigma(n,:).^4;
% 
% k1 = sqrt(term1+sqrt(term2));
% k2 = sqrt(term1-sqrt(term2));
% 
% 
% figure(100)
% plot(h,term1); hold on;
% plot(h,term2); 
% legend('term1','term2','Location','NorthWest'); legend boxoff
% xlabel('\alpha')
% title(['r=',num2str(R(n))])
% 
% figure(200)
% x0=10; y0=10; width=1000; height=400; 
% set(gcf,'position',[x0,y0,width,height])
% 
% subplot(1,2,1)
% plot(h,real(k1)); hold on;
% plot(h,real(k2));
% legend('real(k1)','real(k2)'); legend boxoff
% xlabel('h')
% title(['r=',num2str(R(n))])
% 
% subplot(1,2,2)
% plot(h,imag(k1)); hold on;
% plot(h,imag(k2));
% legend('imag(k1)','imag(k2)'); legend boxoff
% xlabel('h')
% title(['r=',num2str(R(n))])

% Mode or DRV plot: fixed r

% r=0.1

load('mode_or_drv/mode_or_drv_h1_vs_h2_fixed_r01_N200_L8pi_t200.mat')



% % set NaN for all negative growth rates
% 
% mode_or_drv(sigma<0.09) = NaN;

% if zero-local max set to DRV (lies at the boundary)

for ii = 1:length(h1)
    for jj = 1:length(h2)
ind = islocalmax(ww(ii,jj,:)); 
ind(ww(ii,jj,:)<0) = 0; 

if sum(ind)==0
    mode_or_drv(ii,jj) = 1;
end
    end
end

% define the pv gradients based on h1 and h2

q1y = 1-h1; q2y = -1+h2;

[Q1y,Q2y] = meshgrid(q1y,q2y);


val = [1,0.9,0.8,0.7];
% Plot mode vs DRV
subplot(1,3,2)
contourf(Q2y,Q1y,mode_or_drv,[-1 0 2],'Edgecolor','none'); hold on
contour(Q2y,Q1y,sigma,val,'k--','ShowText','on')
%colormap redblue
title('\rm r=0.1')
%legend('Location','eastoutside')
%set(gca,'Ydir','reverse')
xlabel('q_{2y}'); ylabel('q_{1y}')
%title(['R=',num2str(R)])
set(gca,'FontSize',12)
%caxis([-5 5])
map = [200/255 197/255 252/255;255/255 198/255 199/255;  1  1 1]; colormap (map);

annotation('textbox', [0.36, 0.98, 0, 0], 'string', '(b)','fontsize',12)


%brighten(0.8)

%r=0.01
load('mode_or_drv/mode_or_drv_h1_vs_h2_fixed_r001_N200_L8pi_t200.mat')


% set NaN for all negative growth rates

%mode_or_drv(sigma<0.9) = NaN;

% if zero-local max set to DRV (lies at the boundary)

for ii = 1:length(h1)
    for jj = 1:length(h2)
ind = islocalmax(ww(ii,jj,:)); 
ind(ww(ii,jj,:)<0) = 0; 

if sum(ind)==0
    mode_or_drv(ii,jj) = 1;
end
    end
end

% define the pv gradients based on h1 and h2

q1y = 1-h1; q2y = -1+h2;

[Q1y,Q2y] = meshgrid(q1y,q2y);

val = [0.95,1.1,1.3,1.5];
% Plot mode vs DRV
subplot(1,3,3)
contourf(Q2y,Q1y,mode_or_drv,[-1 0 2],'Edgecolor','none'); hold on
contour(Q2y,Q1y,sigma,val,'k--','ShowText','on')
%colormap redblue
title('\rm r=0.01')
%legend('Location','eastoutside')
%set(gca,'Ydir','reverse')
xlabel('q_{2y}'); ylabel('q_{1y}')
%title(['R=',num2str(R)])
set(gca,'FontSize',12)
%caxis([-5 5])
map = [200/255 197/255 252/255;255/255 198/255 199/255;  1  1 1]; colormap (map);

annotation('textbox', [0.64, 0.98, 0, 0], 'string', '(c)','fontsize',12)

%saveas(gcf,'/net/halo/disk28/disk7/mkohl/Mode_Kohl_OGorman/DRV_Paper/Figures/periodic_wave_or_drv','epsc')





