% plot the beta rossby number as suggested by MM:
% beta_rossby = L^2*qy/U, use different scales for L

close all; clear;

figure(1)

x0=10; y0=10; width= 7/7*1100; height=7/7*300; 
set(gcf,'position',[x0,y0,width,height])

set(gcf,'position',[x0,y0,width,height])


load('mode_or_drv/mode_or_drv_N200_L8pi_t200_with_stream.mat')
Q = zeros(size(ww));
Qx = zeros(size(ww));
V = zeros(size(ww));


dx = L/N;
qy = 1-h;

[Qy,RR] = meshgrid(qy,R);
ratio = zeros(size(Qy));

for ii = 1:size(ww,1)
    for jj = 1:size(ww,2)
 phi = squeeze(phi_stream(ii,jj,:));
 tau = squeeze(tau_stream(ii,jj,:));
 w = squeeze(ww(ii,jj,:));
 q = d_2x(N,dx)*(phi-tau)+tau;
 v = d_1x(N,dx)*(phi-tau);
 qx = d_3x(N,dx)*(phi-tau)+d_1x(N,dx)*tau;
 
 pholder = abs(v*Qy(ii,jj))./abs(qx);
 
 [test,ind] =min(abs(find(w==max(w))-find(w<0)));
 
 
 pb = v*Qy(ii,jj);
 bp = qx;
 
 
 indl = 3;
ratio(ii,jj) = norm(abs(pb(w<0)))/norm(abs(bp(w<0)));
%  if ind<=indl
%     ratio(ii,jj) = NaN;
%  else
%  ratio(ii,jj) = rms(abs(pb(ind-indl:ind))./abs(bp(ind-indl:ind)));
%  end
%  ratio(ii,jj) = rms(pholder(w<0));
%  ratio(ii,jj) = rms(pholder(w<0));
 
%  Q(ii,jj,:) = d_2x(N,dx)*(phi-tau)+tau;
%  Qx(ii,jj,:) = d_3x(N,dx)*(phi-tau)+d_1x(N,dx)*tau;
%  V(ii,jj,:) = d_1x(N,dx)*(phi-tau);
    end
end

ratio(sigma<0.09) = NaN;
%ratio = ratio -1;

%ratio = log(ratio);
ratio = ratio -1;

figure(10)
contourf(RR,Qy,ratio,'Edgecolor','none');  
%%contourf(RR,Qy,beta_rossby,'Edgecolor','none');  
%colormap(redblue(15));
colormap(flipud(redblue(40)));
colorbar; 
%caxis([-6 6]); 
%caxis([-max(beta_rossby(:)) max(beta_rossby(:))])
%caxis([min(ratio(:)) -min(ratio(:))])
caxis([-4 4])
xlabel('r'); ylabel('q_{1y}=-q_{2y}')
set(gca,'FontSize',12)
title('\rm R_{\beta} - 1')
annotation('textbox', [0.62, 0.77, 0.2, 0], 'string', 'Periodic Wave','Edgecolor','none','fontsize',12)
annotation('textbox', [0.2, 0.55, 0, 0], 'string', 'DRV','fontsize',12)
annotation('textbox', [0.55, 0.35, 0, 0], 'string', 'Stable','fontsize',12)

%saveas(gcf,'/net/halo/disk28/disk7/mkohl/Mode_Kohl_OGorman_with_drag/figures_reviewers/beta_rossby_number','epsc')

% set all the stable modes to NaN values

mode_or_drv(sigma<0.09) = 2;

% if zero-local max set to DRV (lies at the boundary)

% for ii = 1:length(R)
%     for jj = 1:length(h)
% ind = islocalmax(ww(ii,jj,:)); 
% ind(ww(ii,jj,:)<0) = 0; 
% 
% if sum(ind)==0
%     mode_or_drv(ii,jj) = 1;
% end
%     end
% end

% dx = L/N;
% for ii = 1:length(R)
%     for jj = 1:length(h)
%         w = squeeze(ww(ii,jj,:));
%         if mode_or_drv(ii,jj) ==1
%             LL(ii,jj) = 1/sigma(ii,jj);
%         elseif mode_or_drv(ii,jj) == -1
%             dist = abs(find(w==min(w))-find(w>0));
%             dist = sort(dist);
%             LL(ii,jj) = 2*dx*dist(1);
%         else
%             LL(ii,jj) = NaN;
%         end
%     end
% end
k = [0:N/2,-N/2+1:-1]*2*pi/L;
indx = N/2;
dx = L/N;

for ii = 1:size(ww,1)
    for jj = 1:size(ww,2)
      phi = squeeze(phi_stream(ii,jj,:));
      tau = squeeze(tau_stream(ii,jj,:));
      v2 = d_1x(N,dx)*(phi-tau);
      q = d_2x(N,dx)*(phi-tau)+tau;
      w = squeeze(ww(ii,jj,:));
      power = abs(fft(q)).^2;
      %[pholder,scale] = max(power);
      scale = sum(k(1:indx)'.*power(1:indx))/sum(power(1:indx));
      LL(ii,jj) = 2*pi/scale;
    end
end

LL = 1./sigma;
LL(sigma<0.09) = NaN;

      
qy = 1-h;

[Qy,RR] = meshgrid(qy,R);

QQy = repmat(Qy,1,1,size(ww,3));

test = abs(V.*QQy)./abs(Qx);

L_descent = (-sigma-sqrt(sigma.^2-4*abs(Qy)))./(2*abs(Qy));

levels = [1.4,1.2,0.8,0.6,0.4,0.09];

%L = 1./sigma;
beta_rossby = LL.^2.*abs(Qy)./1;
%beta_rossby = LL.^2.*abs(Qy)./1-1;
beta_rossby(sigma<0.09) = NaN;

beta_rossby = beta_rossby -1;
%beta_rossby = log(beta_rossby);
L(sigma<0.09) = NaN;

figure(2)

contourf(RR,Qy,beta_rossby,'Edgecolor','none');  
colormap(redblue(20));
colorbar; 
%caxis([-6 6]); 
%caxis([-max(beta_rossby(:)) max(beta_rossby(:))])
caxis([min(beta_rossby(:)) -min(beta_rossby(:))])
xlabel('r'); ylabel('q_{1y}=-q_{2y}')
set(gca,'FontSize',12)
title('R_{\beta} - 1')

clear('LL')










subplot(1,3,1)

contourf(RR,Qy,log(beta_rossby),'Edgecolor','none');  
%%contourf(RR,Qy,beta_rossby,'Edgecolor','none');  
%colormap(redblue(15));
colormap redblue;
colorbar; 
caxis([-6 6]); 
%caxis([-max(beta_rossby(:)) max(beta_rossby(:))])
caxis([min(beta_rossby(:)) -min(beta_rossby(:))])
xlabel('r'); ylabel('q_{1y}')
set(gca,'FontSize',12)
title('\rm q_{1y} = -q_{2y}')

clear('LL')

% Mode or DRV plot: fixed r

% r=0.1

load('mode_or_drv/mode_or_drv_h1_vs_h2_fixed_r01_N200_L8pi_t200.mat')

% if zero-local max set to DRV (lies at the boundary)

for ii = 1:size(ww,1)
    for jj = 1:size(ww,2)
ind = islocalmax(ww(ii,jj,:)); 
ind(ww(ii,jj,:)<0) = 0; 

if sum(ind)==0
    mode_or_drv(ii,jj) = 1;
end
    end
end

for ii = 1:size(ww,1)
    for jj = 1:size(ww,2)
        w = squeeze(ww(ii,jj,:));
        if mode_or_drv(ii,jj) ==1
            LL(ii,jj) = 1/sigma(ii,jj);
        elseif mode_or_drv(ii,jj) == -1
            dist = abs(find(w==min(w))-find(w>0));
            dist = sort(dist);
            LL(ii,jj) = 2*dx*dist(1);
        else
            LL(ii,jj) = NaN;
        end
    end
end

% k = [0:N/2,-N/2+1:-1]*2*pi/L;
% indx = N/2;
% for ii = 1:size(ww,1)
%     for jj = 1:size(ww,2)
%         w = squeeze(ww(ii,jj,:));
%       power = abs(fft(w)).^2;
%       scale = sum(k(1:indx)'.*power(1:indx))/sum(power(1:indx));
%       LL(ii,jj) = 2*pi/scale;
%     end
% end
% LL(sigma<0.09) = NaN;



% define the pv gradients based on h1 and h2

q1y = 1-h1; q2y = -1+h2;

[Q1y,Q2y] = meshgrid(q1y,q2y);

for ii = 1:size(Q1y,1)
    for jj = 1:size(Q1y,2)
        Qymax(ii,jj) = max(abs(Q1y(ii,jj)),abs(Q2y(ii,jj)));
    end
end


val = [1,0.9,0.8,0.7];

%L = 1./sigma;
beta_rossby = LL.^2.*abs(Qymax)./1;
%beta_rossby = LL.^2.*abs(Qymax)./1-1;
beta_rossby(sigma<0.09) = NaN;

subplot(1,3,2)

contourf(Q2y,Q1y,log(beta_rossby),'Edgecolor','none');  
%contourf(Q2y,Q1y,beta_rossby,'Edgecolor','none');  
%colormap(redblue(15));
colormap redblue;
colorbar; 
caxis([-6 6]); 
%caxis([-max(beta_rossby(:)) max(beta_rossby(:))])
caxis([min(beta_rossby(:)) -min(beta_rossby(:))])
xlabel('q_{2y}'); ylabel('q_{1y}')
set(gca,'FontSize',12)
title('\rm r=0.1')

clear('LL');

%r=0.01
load('mode_or_drv/mode_or_drv_h1_vs_h2_fixed_r001_N200_L8pi_t200.mat')

% if zero-local max set to DRV (lies at the boundary)

for ii = 1:size(ww,1)
    for jj = 1:size(ww,2)
ind = islocalmax(ww(ii,jj,:)); 
ind(ww(ii,jj,:)<0) = 0; 

if sum(ind)==0
    mode_or_drv(ii,jj) = 1;
end
    end
end

for ii = 1:size(ww,1)
    for jj = 1:size(ww,2)
        w = squeeze(ww(ii,jj,:));
        
        if mode_or_drv(ii,jj) ==1
            LL(ii,jj) = 1/sigma(ii,jj);
        elseif mode_or_drv(ii,jj) == -1
            dist = abs(find(w==min(w))-find(w>0));
            dist = sort(dist);
            LL(ii,jj) = 2*dx*dist(1);
        else
            LL(ii,jj) = NaN;
        end
    end
end

% k = [0:N/2,-N/2+1:-1]*2*pi/L;
% indx = N/2;
% for ii = 1:size(ww,1)
%     for jj = 1:size(ww,2)
%         w = squeeze(ww(ii,jj,:));
%       power = abs(fft(w)).^2;
%       scale = sum(k(1:indx)'.*power(1:indx))/sum(power(1:indx));
%       LL(ii,jj) = 2*pi/scale;
%     end
% end
% LL(sigma<0.09) = NaN;

% define the pv gradients based on h1 and h2

q1y = 1-h1; q2y = -1+h2;

[Q1y,Q2y] = meshgrid(q1y,q2y);

for ii = 1:size(Q1y,1)
    for jj = 1:size(Q1y,2)
        Qymax(ii,jj) = max(abs(Q1y(ii,jj)),abs(Q2y(ii,jj)));
    end
end


val = [1,0.9,0.8,0.7];

%L = 1./sigma;
beta_rossby = LL.^2.*abs(Qymax)./1;
%beta_rossby = LL.^2.*abs(Qymax)./1-1;
beta_rossby(sigma<0.09) = NaN;

subplot(1,3,3)

contourf(Q2y,Q1y,log(beta_rossby),'Edgecolor','none');  
%contourf(Q2y,Q1y,beta_rossby,'Edgecolor','none');  
%colormap(redblue(15));
caxis([-6 6]); 
%caxis([-max(beta_rossby(:)) max(beta_rossby(:))])
caxis([min(beta_rossby(:)) -min(beta_rossby(:))])
colormap redblue;
colorbar; 

xlabel('q_{2y}'); ylabel('q_{1y}')
set(gca,'FontSize',12)
title('\rm r=0.01')

annotation('textbox', [0.64, 0.98, 0, 0], 'string', '(c)','fontsize',12)




