% Plot the PV Budget for the DRV example for geostrophic vs.
% semigeostrophic flow

close all; clear;

load('DRV_Paper/DRV_example.mat')

%% Semi-geostrophic Theory

epsilon = 1; % Rossby-Number

res = 120;%140;
phi_end = phi_end/res; tau_end = tau_end/res;

% Define the meridional velocity in top (1) and bottom (2) layer

v1 = d_1x(N,dx)*(phi_end+tau_end); v1 = [v1;v1(1)];
v2 = d_1x(N,dx)*(phi_end-tau_end); v2 = [v2;v2(1)];

Zeta1 = d_2x(N,dx)*(phi_end+tau_end);
Zeta2 = d_2x(N,dx)*(phi_end-tau_end);

% Apply the coordinate transformation x_i = X-v_i to plot the variables
% in physical space

x1 = x-epsilon*v1'; 
x2 = x-epsilon*v2'; 

phi1 = (phi_end+tau_end)-epsilon*1/2*v1(1:end-1).^2;
phi2 = (phi_end-tau_end)-epsilon*1/2*v2(1:end-1).^2;

% test: should be the same Phi_X = phi_x

v1s = d_1xs(N,x1)*phi1; v1s = [v1s;v1s(1)];
v2s = d_1xs(N,x2)*phi2; v2s = [v2s;v2s(1)];

% Now transform the terms of the PV(X) = Phi_iXX +- Tau(X)
% Phi_iXX = phi_ixx/(1+phi_ixx)
% Tau(X) = tau(x+v)

phi1xx = d_2xs(N,x1)*phi1; phi2xx = d_2xs(N,x2)*phi2;
% 
% zeta1 = d_1xs(N,x1)*v1(1:end-1); 
% 
% zeta2 = d_1xs(N,x2)*v2(1:end-1);

J1 = ones(N,1)+d_2xs(N,x1)*phi1;
J1m = ones(N,1)-d_2x(N,dx)*(phi_end+tau_end);

J2 = ones(N,1)+d_2xs(N,x2)*phi2;
J2m = ones(N,1)-d_2x(N,dx)*(phi_end-tau_end);


PV1X = Zeta1-tau_end; PV2X = Zeta2+tau_end;

PV1 = J1.*PV1X; PV1 = [PV1;PV1(1)];
PV2 = J2.*PV2X; PV2 = [PV2;PV2(2)];


PV1X = [PV1X;PV1X(1)];
PV2X = [PV2X;PV2X(1)];

RHS = 2*d_1x(N,dx)*Phi(:,end)/res-d_1x(N,dx)*phi_end;
% 
w = Omega_Solver(RHS,R,N,dx);

r = r_factor(w,R);

PV1XG = -(1-r).*w; PV1XG = [PV1XG;PV1XG(1)];

PV2XG = (1-r).*w; PV2XG = [PV2XG;PV2XG(1)];

PV1G = -J1.*(1-r).*w; PV1G = [PV1G;PV1G(1)];

PV2G = J2.*(1-r).*w; PV2G = [PV2G;PV2G(1)];


AmpX = max(abs(PV2X))/max(abs(PV1X));
Amp = max(abs(PV2))/max(abs(PV1));

% Plot X-space

hfig = figure(1);
pos = get(hfig,'position');
set(hfig,'position',pos.*[.5 1 2 1])

[val,int] = max(PV2X); intL = int-100; intR = int+100;

subplot(1,2,1)
plot(x(intL:intR),PV2X(intL:intR),'Color','r','linewidth',2.5); hold on;
plot(x(intL:intR),PV2XG(intL:intR)/2,'Color','r','linewidth',2,'linestyle','--'); hold on;
plot(x(intL:intR),PV1X(intL:intR)+8*max(abs(PV2X)),'Color','b','linewidth',2.5); hold on
plot(x(intL:intR),PV1XG(intL:intR)/2+8*max(abs(PV1X)),'Color','b','linewidth',2,'Linestyle','--');
xlabel('X')
set(gca,'fontsize', 14);
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
xlim([x(intL) x(intR)])
ylim([-0.05 8])
title('Quasigeostrophic')

[val,int] = max(PV2); intL = int-50; intR = int+50;

subplot(1,2,2)
plot(x1(intL:intR),PV1(intL:intR)+8*max(abs(PV2X)),'Color','b','linewidth',2.5); hold on
plot(x1(intL:intR),PV1G(intL:intR)/2+8*max(abs(PV2X)),'Color','b','linewidth',2,'Linestyle','--'); hold on;
plot(x2(intL:intR),PV2(intL:intR),'Color','r','linewidth',2.5); hold on;
plot(x2(intL:intR),PV2G(intL:intR)/2,'Color','r','linewidth',2,'linestyle','--'); 
xlabel('x')
set(gca,'fontsize', 14);
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
xlim([x2(intL) x2(intR)])
ylim([-0.05 8])
title('Semigeostrophic')
% legend('PV_{up}','Diabatic_{up}','PV_{down}','Diabatic_{down}');
% legend boxoff 
% 
% % add a bit space to the figure
% fig = gcf;
% fig.Position(3) = fig.Position(3) + 250;
% % add legend
% Lgnd = legend('show');
% Lgnd.Position(1) = 0.01;
% Lgnd.Position(2) = 0.4;
%saveas(gcf,'DRV_Paper/Figures/DRV_qg_vs_sg','epsc')

% figure(3)
% plot(x2,PV2/max(PV2),'linewidth',2.5); hold on;
% plot(x2,v2,'linewidth',2.5); hold on;
% %plot(x2,[w;w(1)]/max(w),'linewidth',2.5)
% legend('PV2','v2'); legend boxoff
% xlabel('x')
% xlim([x2(1) x2(end)])
% title('Bottom Layer')
% 
% figure(4)
% plot(x1,PV1/max(PV1),'linewidth',2.5); hold on;
% plot(x1,v1,'linewidth',2.5); hold on;
% %plot(x2,[w;w(1)]/max(w),'linewidth',2.5)
% legend('PV1','v1'); legend boxoff
% xlabel('x')
% xlim([x1(1) x1(end)])
% title('Top Layer')
% 
% 
% [val,ind] = max(abs(x-x2)); indL = ind-20; indR = ind+20;
% 
% figure(5)
% plot(x(indL:indR),x(indL:indR),'linewidth',1.5); hold on;
% plot(x(indL:indR),x1(indL:indR),'linewidth',1.5); hold on;
% plot(x(indL:indR),x2(indL:indR),'linewidth',1.5); 
% xlabel('X')
% legend('X','x1','x2','Location','Northwest'); legend boxoff
% xlim([x(indL) x(indR)])
% title('X-x Transform')


% Plot Hoskins75 equation (21) PV = -(Phi1-Phi2)/(1-Phi1_XX)
Phi1 = phi_end+tau_end; 
Phi2 = phi_end-tau_end;

Q1 = (1-(Phi1-Phi2))./(1-Zeta1); Q1 = [Q1;Q1(1)]-1;
Q2 = (1+(Phi1-Phi2))./(1-Zeta2); Q2 = [Q2;Q2(1)]-1;

q1 = phi1xx - (phi1-phi2); q1 = [q1;q1(1)];
q2 = phi2xx + (phi1-phi2); q2 = [q2;q2(1)];

figure(2)
plot(x,Q2,'Color','r'); hold on; plot(x,Q1+3/2*max(Q2),'Color','b'); hold on;
plot(x,q2,'Color','r','Linestyle','--'); hold on; plot(x,q1+3/2*max(Q2),'Color','b','Linestyle','--');
xlim([x(1) x(end)])
legend('Ertel','','Pseudo',''); legend boxoff
