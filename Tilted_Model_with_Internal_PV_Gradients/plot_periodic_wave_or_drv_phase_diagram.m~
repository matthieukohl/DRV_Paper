% Mode or DRV plot: Ascent lenghts

close all; clear;

x0=10; y0=10; width= 8/7*1100; height=2*8/7*300; 
set(gcf,'position',[x0,y0,width,height])

load('mode_or_drv/mode_or_drv_N200_L8pi_t200_new.mat')
%load('mode_or_drv/mode_or_drv_N50_L8pi_t200.mat')

% set all the stable modes to NaN values

mode_or_drv(sigma<0.09) = 2;
%b(sigma<0.09) = NaN;
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

[RR,hh] = meshgrid(qy,R);

levels = linspace(min(b(:)),max(b(:)),8);
levels = round(levels,1);

% calculate b

for ii = 1:size(ww,1)
    for jj = 1:size(ww,2)
  w = squeeze(ww(ii,jj,:));
diff = abs(find(w==max(w))-find(w<0));
diff = sort(diff);
dx = L/N;
b(ii,jj) = 0.5*(dx*(diff(1)+diff(2)));
    end
end

b(sigma<0.09) = NaN;

levels = [1.4,1.2,0.8,0.6,0.4,0.09];
s1 = subplot(2,3,1);


contourf(hh,RR,mode_or_drv,[-1 0 2],'EdgeColor','none'); hold on;
map = [200/255 197/255 252/255;255/255 198/255 199/255;  1  1 1];% 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1]; 
colormap(s1,map);
caxis([-1 2]); % this is needed so that three colors only span variable mode_or_drv
contour(hh,RR,sigma,levels,'k--','ShowText','on'); hold on;
t = title('\rm q_{2y} = -q_{1y}');
set(t,'position',get(t,'position')+[0 0.15 0])
%xlabel('r'); 
ylabel('q_{1y}')
set(gca,'FontSize',12)

annotation('textbox', [0.07, 0.96, 0, 0], 'string', '(a)','fontsize',12)
annotation('textbox', [0.20, 0.92, 0.2, 0], 'string', 'Periodic Wave','Edgecolor','none','fontsize',12)
annotation('textbox', [0.24, 0.83, 0, 0], 'string', 'DRV','fontsize',12)
annotation('textbox', [0.27, 0.67, 0, 0], 'string', 'Stable','fontsize',12)



s4 = subplot(2,3,4);
dispy = 0.03;
set(s4,'position',get(s4,'position')+[0 dispy 0 0]);

levels = linspace(min(b(:))+0.2,max(b(:)),6);
levels = round(levels,1);

contourf(hh,RR,mode_or_drv,[-1 0 2],'EdgeColor','none'); hold on;
map = [200/255 197/255 252/255;255/255 198/255 199/255;  1  1 1];% 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1]; 
colormap(s4,map);
caxis([-1 2]); % this is needed so that three colors only span variable mode_or_drv
contour(hh,RR,b,levels,'k--','ShowText','on'); hold on;
%title('\rm q_{2y} = -q_{1y}')
xl = xlabel('r'); 
ylabel('q_{1y}')
set(gca,'FontSize',12)
set(xl,'position',get(xl,'position')+[0 -0.15 0])

annotation('textbox', [0.07, 0.51, 0, 0], 'string', '(d)','fontsize',12)


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

% calculate b

for ii = 1:size(ww,1)
    for jj = 1:size(ww,2)
  w = squeeze(ww(ii,jj,:));
diff = abs(find(w==max(w))-find(w<0));
diff = sort(diff);
dx = L/N;
b(ii,jj) = 0.5*(dx*(diff(1)+diff(2)));
    end
end

% define the pv gradients based on h1 and h2

q1y = 1-h1; q2y = -1+h2;

[Q1y,Q2y] = meshgrid(q1y,q2y);

val = [1,0.9,0.8,0.7];

s2 = subplot(2,3,2);
contourf(Q2y,Q1y,mode_or_drv,[-1 0 2],'Edgecolor','none'); hold on
contour(Q2y,Q1y,sigma,val,'k--','ShowText','on')
t = title('\rm r=0.1');
set(t,'position',get(t,'position')+[0 0.15 0])
%xlabel('q_{2y}'); ylabel('q_{1y}')
set(gca,'FontSize',12)
map = [200/255 197/255 252/255;255/255 198/255 199/255;  1  1 1]; 
colormap (s2,map);
caxis([-1 2])

annotation('textbox', [0.35, 0.96, 0, 0], 'string', '(b)','fontsize',12)
xticks([-1 -0.5 0 0.5 1])
xticklabels({'-1','-0.5','0','0.5','1'})


val = linspace(min(b(:)),max(b(:)),7);
val = round(val,2);

% Plot mode vs DRV
s5 = subplot(2,3,5);
set(s5,'position',get(s5,'position')+[0 dispy 0 0])
contourf(Q2y,Q1y,mode_or_drv,[-1 0 2],'Edgecolor','none'); hold on
contour(Q2y,Q1y,b,val,'k--','ShowText','on')
%title('\rm r=0.1')
xl = xlabel('q_{2y}'); %ylabel('q_{1y}')
set(xl,'position',get(xl,'position')+[0 -0.15 0])
set(gca,'FontSize',12)
map = [200/255 197/255 252/255;255/255 198/255 199/255;  1  1 1]; 
colormap (s5,map);
caxis([-1 2])
annotation('textbox', [0.35, 0.51, 0, 0], 'string', '(e)','fontsize',12)
xticks([-1 -0.5 0 0.5 1])
xticklabels({'-1','-0.5','0','0.5','1'})

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

% calculate b

for ii = 1:size(ww,1)
    for jj = 1:size(ww,2)
  w = squeeze(ww(ii,jj,:));
diff = abs(find(w==max(w))-find(w<0));
diff = sort(diff);
dx = L/N;
b(ii,jj) = 0.5*(dx*(diff(1)+diff(2)));
    end
end

% define the pv gradients based on h1 and h2

q1y = 1-h1; q2y = -1+h2;

[Q1y,Q2y] = meshgrid(q1y,q2y);


val = [0.95,1.1,1.3,1.5];

% Plot mode vs DRV
s3 = subplot(2,3,3);

contourf(Q2y,Q1y,mode_or_drv,[-1 1],'Edgecolor','none'); hold on
contour(Q2y,Q1y,sigma,val,'k--','ShowText','on')
t = title('\rm r=0.01');
set(t,'position',get(t,'position')+[0 0.15 0])
%xlabel('q_{2y}'); ylabel('q_{1y}')
set(gca,'FontSize',12)
map = [200/255 197/255 252/255;255/255 198/255 199/255;255/255 198/255 199/255;  1  1 1]; 
colormap (s3,map);
caxis([-1 2])
xticks([-1 -0.5 0 0.5 1])
xticklabels({'-1','-0.5','0','0.5','1'})
annotation('textbox', [0.63, 0.96, 0, 0], 'string', '(c)','fontsize',12)


% Plot mode vs DRV
%val = [0.18,0.25,0.4];
val = [0.1885, 0.22];
s6 = subplot(2,3,6);
set(s6,'position',get(s6,'position')+[0 dispy 0 0])

contourf(Q2y,Q1y,mode_or_drv,[-1 1],'Edgecolor','none'); hold on
contour(Q2y,Q1y,b,val,'k--','ShowText','on')
%title('\rm r=0.01')
xl = xlabel('q_{2y}'); %ylabel('q_{1y}')
set(xl,'position',get(xl,'position')+[0 -0.15 0])
set(gca,'FontSize',12)
map = [200/255 197/255 252/255;255/255 198/255 199/255;255/255 198/255 199/255;  1  1 1]; 
colormap (s6,map);
caxis([-1 2])
xticks([-1 -0.5 0 0.5 1])
xticklabels({'-1','-0.5','0','0.5','1'})

annotation('textbox', [0.63, 0.51, 0, 0], 'string', '(f)','fontsize',12)



annotation('textbox', [0.93, 0.77, 0, 0], 'string', 'Growth Rate (\sigma)','fontsize',12)
annotation('textbox', [0.93, 0.37, 0, 0], 'string', 'Half Ascent Length (b)','fontsize',12)

annotation('textbox', [0.483, 0.89, 0.2, 0], 'string', 'Periodic Wave','Edgecolor','none','fontsize',12)
annotation('textbox', [0.54, 0.70, 0, 0], 'string', 'DRV','fontsize',12)

annotation('textbox', [0.77, 0.89, 0.2, 0], 'string', 'Periodic Wave','Edgecolor','none','fontsize',12)
annotation('textbox', [0.80, 0.70, 0, 0], 'string', 'DRV','fontsize',12)


%saveas(gcf,'/net/halo/disk28/disk7/mkohl/Mode_Kohl_OGorman/DRV_Paper/Figures/periodic_wave_or_drv_with_sigma_and_b_version1','epsc')
