% plot the beta rossby number as suggested by MM:
% beta_rossby = L^2*qy/U, use different scales for L

close all; clear;

figure(1)
x0=10; y0=10; width= 8/7*1100; height=8/7*300; 
set(gcf,'position',[x0,y0,width,height])

set(gcf,'position',[x0,y0,width,height])


load('mode_or_drv/mode_or_drv_N200_L8pi_t200_new.mat')

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

[Qy,RR] = meshgrid(qy,R);

levels = [1.4,1.2,0.8,0.6,0.4,0.09];

L = 2*b;
beta_rossby = L.^2.*abs(Qy)./1;
beta_rossby(sigma<0.09) = NaN;
% beta_rossby(beta_rossby<1) = 1;
% beta_rossby(beta_rossby>1) = -1;

L(sigma<0.09) = NaN;

figure(1)
subplot(1,3,1)
% contourf(RR,Qy,beta_rossby,[-1 1],'EdgeColor','None'); hold on;
% colormap redblue
% caxis([-5 5])
% %map = [200/255 197/255 252/255;255/255 198/255 199/255;  1  1 1]; colormap (map);
% contour(RR,Qy,L,'k','Linewidth',1.5,'ShowText','on')
% title('\rm L = 2b')
% %set(gca,'Ydir','reverse')
% xlabel('r'); %ylabel('q_{1y}=-q_{2y}')
% ylabel('q_{1y}')
% set(gca,'FontSize',12)

contourf(RR,Qy,log(beta_rossby),'Edgecolor','none');  
colormap(redblue(15));
colorbar; 
caxis([-6 6]); 
xlabel('r'); ylabel('q_{1y}')
set(gca,'FontSize',12)
title('\rm L = 2b')

L = 1./sigma;
beta_rossby = L.^2.*abs(Qy)./1;
beta_rossby(sigma<0.09) = NaN;
%beta_rossby(beta_rossby<1) = 1;
%beta_rossby(beta_rossby>1) = -1;

L(sigma<0.09) = NaN;

subplot(1,3,2)
% contourf(RR,Qy,beta_rossby,[-1 1],'EdgeColor','None'); hold on;
% colormap redblue
% caxis([-5 5])
% %map = [200/255 197/255 252/255;255/255 198/255 199/255;  1  1 1]; colormap (map);
% contour(RR,Qy,L,'k','Linewidth',1.5,'ShowText','on')
% title('\rm L = 1/\sigma')
% %set(gca,'Ydir','reverse')
% xlabel('r'); %ylabel('q_{1y}=-q_{2y}')
% ylabel('q_{1y}')
% set(gca,'FontSize',12)

contourf(RR,Qy,log(beta_rossby),'Edgecolor','none');  
colormap(redblue(15));
colorbar; 
caxis([-6 6]); 
xlabel('r'); ylabel('q_{1y}')
set(gca,'FontSize',12)
title('\rm L = 1/\sigma')

L = 2*b+1./sigma;
beta_rossby = L.^2.*abs(qy)./1;
beta_rossby(sigma<0.09) = NaN;
%beta_rossby(beta_rossby<1) = 1;
%beta_rossby(beta_rossby>1) = -1;

L(sigma<0.09) = NaN;

subplot(1,3,3)
% contourf(RR,Qy,beta_rossby,[-1 1],'EdgeColor','None'); hold on;
% colormap redblue
% caxis([-5 5])
% %map = [200/255 197/255 252/255;255/255 198/255 199/255;  1  1 1]; colormap (map);
% contour(RR,Qy,L,'k','Linewidth',1.5,'ShowText','on')
% title('\rm L = 1/\sigma+2b')
% %set(gca,'Ydir','reverse')
% xlabel('r'); %ylabel('q_{1y}=-q_{2y}')
% ylabel('q_{1y}')
% set(gca,'FontSize',12)

contourf(RR,Qy,log(beta_rossby),'Edgecolor','none');  
colormap(redblue(15));
colorbar;
caxis([-6 6]); 


xlabel('r'); ylabel('q_{1y}')
set(gca,'FontSize',12)
title('\rm L = 1/\sigma+2b')

figure(2)

x0=10; y0=10; width= 8/7*1100; height=8/7*300; 
set(gcf,'position',[x0,y0,width,height])

set(gcf,'position',[x0,y0,width,height])


load('mode_or_drv/mode_or_drv_N200_L8pi_t200_new.mat')

% set all the stable modes to NaN values

mode_or_drv(sigma<0.09) = 2;

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

[Qy,RR] = meshgrid(qy,R);

levels = [1.4,1.2,0.8,0.6,0.4,0.09];

L = 1./sigma;
beta_rossby = L.^2.*abs(Qy)./1;
beta_rossby(sigma<0.09) = NaN;


L(sigma<0.09) = NaN;

subplot(1,3,1)

contourf(RR,Qy,log(beta_rossby),'Edgecolor','none');  
colormap(redblue(15));
colorbar; 
caxis([-6 6]); 
xlabel('r'); ylabel('q_{1y}')
set(gca,'FontSize',12)
title('\rm q_{1y} = -q_{2y}')

% Mode or DRV plot: fixed r

% r=0.1

load('mode_or_drv/mode_or_drv_h1_vs_h2_fixed_r01_N200_L8pi_t200.mat')

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

for ii = 1:size(Q1y,1)
    for jj = 1:size(Q1y,2)
        Qymax(ii,jj) = max(abs(Q1y(ii,jj)),abs(Q2y(ii,jj)));
    end
end


val = [1,0.9,0.8,0.7];

L = 1./sigma;
beta_rossby = L.^2.*abs(Qymax)./1;
beta_rossby(sigma<0.09) = NaN;

subplot(1,3,2)

contourf(Q2y,Q1y,log(beta_rossby),'Edgecolor','none');  
colormap(redblue(15));
colorbar; 
caxis([-6 6]); 
xlabel('q_{2y}'); ylabel('q_{1y}')
set(gca,'FontSize',12)
title('\rm r=0.1')


%r=0.01
load('mode_or_drv/mode_or_drv_h1_vs_h2_fixed_r001_N200_L8pi_t200.mat')

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

for ii = 1:size(Q1y,1)
    for jj = 1:size(Q1y,2)
        Qymax(ii,jj) = max(abs(Q1y(ii,jj)),abs(Q2y(ii,jj)));
    end
end


val = [1,0.9,0.8,0.7];

L = 1./sigma;
beta_rossby = L.^2.*abs(Qymax)./1;
beta_rossby(sigma<0.09) = NaN;

subplot(1,3,3)

contourf(Q2y,Q1y,log(beta_rossby),'Edgecolor','none');  
colormap(redblue(15));
colorbar; 
caxis([-6 6]); 
xlabel('q_{2y}'); ylabel('q_{1y}')
set(gca,'FontSize',12)
title('\rm r=0.01')

annotation('textbox', [0.64, 0.98, 0, 0], 'string', '(c)','fontsize',12)




