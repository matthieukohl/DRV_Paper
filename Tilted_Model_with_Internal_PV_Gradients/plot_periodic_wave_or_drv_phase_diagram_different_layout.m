% Mode or DRV plot: include extra row for growth rates and ascent lenghts

close all; clear;

x0=10; y0=10; width= 8/7*1100; height=2*8/7*300; 
set(gcf,'position',[x0,y0,width,height])
margins1 = [0.14,0.10];
margins2 = margins1;

load('mode_or_drv/mode_or_drv_N200_L8pi_t200_new.mat')
%load('mode_or_drv/mode_or_drv_N200_L8pi_t200_drag0.33.mat')
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

%s1 = subplot(2,3,1);
s1 = subplot_tight(2,3,1,margins1);

contourf(hh,RR,mode_or_drv,[-1 0 2],'EdgeColor','none'); 
%contour(hh,RR,sigma,levels,'k--','ShowText','on'); hold on; 
%contour(hh,RR,b,'g:','ShowText','on');
map = [200/255 197/255 252/255;255/255 198/255 199/255;  1  1 1]; 
colormap(s1,map);
cmp = colormap;
t = title('\rm q_{2y} = -q_{1y}');
set(t,'position',get(t,'position')+[0 0.15 0])
xlabel('r'); %ylabel('q_{1y}=-q_{2y}')
ylabel('q_{1y}')
set(gca,'FontSize',12)
annotation('textbox', [0.02, 0.89, 0, 0], 'string', '(a)','fontsize',12)
annotation('textbox', [0.14, 0.87, 0.2, 0], 'string', 'Periodic Wave','Edgecolor','none','fontsize',12)
annotation('textbox', [0.14, 0.73, 0, 0], 'string', 'DRV','fontsize',12)
annotation('textbox', [0.23, 0.65, 0, 0], 'string', 'Stable','fontsize',12)
set(s1,'position',get(s1,'position')+[0 0.01 0 0]);




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
sigma(sigma<0.09) = NaN;

%s4 = subplot(2,3,4);
s4 = subplot_tight(2,3,4,margins2);
contourf(hh,RR,sigma,'EdgeColor','none'); hold on; 
contour(hh,RR,b,'k--','ShowText','on');
colormap(s4,autumn(10));
%colorbar
%title('\rm q_{2y} = -q_{1y}')
xlabel('r'); %ylabel('q_{1y}=-q_{2y}')
ylabel('q_{1y}')
set(gca,'FontSize',12)
annotation('textbox', [0.02, 0.49, 0, 0], 'string', '(d)','fontsize',12)
caxis([0.09 1.5])
set(s4,'position',get(s4,'position')+[0 0.04 0 0]);


% Mode or DRV plot: fixed r

% r=0.1

load('mode_or_drv/mode_or_drv_h1_vs_h2_fixed_r01_N200_L8pi_t200.mat')
%load('mode_or_drv/mode_or_drv_h1_vs_h2_fixed_r01_N200_L8pi_t200_drag0.33.mat')



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


% Plot mode vs DRV
%s2 = subplot(2,3,2);
s2 = subplot_tight(2,3,2,margins1);
contourf(Q2y,Q1y,mode_or_drv,[-1 0 2],'Edgecolor','none'); hold on
%colormap redblue
t = title('\rm r=0.1');
set(t,'position',get(t,'position')+[0 0.15 0])
%legend('Location','eastoutside')
%set(gca,'Ydir','reverse')
xlabel('q_{2y}'); ylabel('q_{1y}')
%title(['R=',num2str(R)])
set(gca,'FontSize',12)
%caxis([-5 5])
map = [200/255 197/255 252/255;255/255 198/255 199/255;  1  1 1]; 
colormap(s2,map);
caxis([-1 2])
set(s2,'position',get(s2,'position')+[0 0.01 0 0]);

annotation('textbox', [0.32, 0.89, 0, 0], 'string', '(b)','fontsize',12)

xticks([-1 -0.5 0 0.5 1])
xticklabels({'-1','-0.5','0','0.5','1'})

for ii = 1:size(ww,1)
    for jj = 1:size(ww,2)
  w = squeeze(ww(ii,jj,:));
diff = abs(find(w==max(w))-find(w<0));
diff = sort(diff);
dx = L/N;
b(ii,jj) = 0.5*(dx*(diff(1)+diff(2)));
    end
end
val = linspace(min(b(:)),max(b(:)),7);
val = round(val,2);

% Plot mode vs DRV
%s5 = subplot(2,3,5);
s5 = subplot_tight(2,3,5,margins2);
contourf(Q2y,Q1y,sigma,'Edgecolor','none'); hold on
contour(Q2y,Q1y,b,val,'k--','ShowText','on');
colormap(s5,autumn(7));
%colorbar
xlabel('q_{2y}'); ylabel('q_{1y}')
set(gca,'FontSize',12)
caxis([0 1.5])
set(s5,'position',get(s5,'position')+[0 0.04 0 0]);

annotation('textbox', [0.32, 0.49, 0, 0], 'string', '(e)','fontsize',12)

xticks([-1 -0.5 0 0.5 1])
xticklabels({'-1','-0.5','0','0.5','1'})



%brighten(0.8)

%r=0.01
load('mode_or_drv/mode_or_drv_h1_vs_h2_fixed_r001_N200_L8pi_t200.mat')
%load('mode_or_drv/mode_or_drv_h1_vs_h2_fixed_r001_N200_L8pi_t200_drag0.33.mat')


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
%s3 = subplot(2,3,3);
s3 = subplot_tight(2,3,3,margins1);
contourf(Q2y,Q1y,mode_or_drv,[-1 0 2],'Edgecolor','none'); hold on
t = title('\rm r=0.01');
set(t,'position',get(t,'position')+[0 0.15 0])
%legend('Location','eastoutside')
%set(gca,'Ydir','reverse')
xlabel('q_{2y}'); ylabel('q_{1y}')
%title(['R=',num2str(R)])
set(gca,'FontSize',12)
%caxis([-5 5])
map = [200/255 197/255 252/255;255/255 198/255 199/255;  1  1 1]; 
caxis([-1 2])
colormap (s3,map);
set(s3,'position',get(s3,'position')+[0 0.01 0 0]);
xticks([-1 -0.5 0 0.5 1])
xticklabels({'-1','-0.5','0','0.5','1'})

annotation('textbox', [0.62, 0.89, 0, 0], 'string', '(c)','fontsize',12)

% Plot mode vs DRV
val = [0.13,0.25];
for ii = 1:size(ww,1)
    for jj = 1:size(ww,2)
  w = squeeze(ww(ii,jj,:));
diff = abs(find(w==max(w))-find(w<0));
diff = sort(diff);
dx = L/N;
b(ii,jj) = 0.5*(dx*(diff(1)+diff(2)));
    end
end


%s6 = subplot(2,3,6);
s6 = subplot_tight(2,3,6,margins2);
contourf(Q2y,Q1y,sigma,'EdgeColor','none'); hold on
contour(Q2y,Q1y,b,val,'k--','ShowText','on')
colormap(s6,autumn(10))

%title('\rm r=0.01')
xlabel('q_{2y}'); ylabel('q_{1y}')
set(gca,'FontSize',12)
annotation('textbox', [0.62, 0.49, 0, 0], 'string', '(f)','fontsize',12)


h = colorbar('south');
caxis([0.09 1.5])
%set(get(h,'title'),'string','\sigma');
%h.Position =[ 0.93    0.1  0.02    0.3];
h.Position =[ 0.3    0.05  0.4    0.02];
set(s6,'position',get(s6,'position')+[0 0.04 0 0]);
xticks([-1 -0.5 0 0.5 1])
xticklabels({'-1','-0.5','0','0.5','1'})

annotation('textbox', [ 0.7    0.06  0.02    0.02], 'string', '\sigma','EdgeColor','none','fontsize',12)


%saveas(gcf,'/net/halo/disk28/disk7/mkohl/Mode_Kohl_OGorman/DRV_Paper/Figures/periodic_wave_or_drv_with_sigma_and_b_version2','epsc')


% growth rate and length scale comparison drag/no-drag


x0=10; y0=10; width= 1000; height= 400; 
set(gcf,'position',[x0,y0,width,height])
subplot(1,2,1)

load('mode_or_drv/mode_or_drv_N200_L8pi_t200_new.mat')
plot(R,sigma(:,7),'b:'); hold on;
plot(R,sigma(:,12),'b'); hold on;
plot(R,sigma(:,15),'b--'); hold on;

load('mode_or_drv/mode_or_drv_N200_L8pi_t200_drag0.33.mat')

plot(R,sigma(:,7),'r:'); hold on;
plot(R,sigma(:,12),'r'); hold on;
plot(R,sigma(:,15),'r--')
xlabel('r');
ylabel('\sigma');

legend('qy=0.5','qy=0','qy=-0.3','with drag','with drag','with drag','Location','NorthEast');
legend boxoff

subplot(1,2,2)

load('mode_or_drv/mode_or_drv_N200_L8pi_t200_new.mat')
plot(R,b(:,7),'b:'); hold on;
plot(R,b(:,12),'b'); hold on;
plot(R,b(:,15),'b--'); hold on;

load('mode_or_drv/mode_or_drv_N200_L8pi_t200_drag0.33.mat')
legend boxoff

plot(R,b(:,7),'r:'); hold on;
plot(R,b(:,12),'r'); hold on;
plot(R,b(:,15),'r--')
xlabel('r');
ylabel('b');

legend('qy=0.5','qy=0','qy=-0.3','with drag','with drag','qy=-0.3 drag','Location','NorthWest');
legend boxoff




