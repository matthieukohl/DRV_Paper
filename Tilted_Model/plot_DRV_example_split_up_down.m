% Plot the DRV-Example for the QJRMS Paper
close all; clear;

load('DRV_Paper/DRV_example.mat')

% Calculate the different terms in the PV Budget

r = r_factor(w,R);

% low layer: sigma*q = Uqx + (1-r)w + mean(r(w)w)

q_low = d_2x(N,dx)*phi_end-d_2x(N,dx)*tau_end+tau_end;

qt_low = sigma*q_low;

uqx_low = (d_3x(N,dx)*phi_end-d_3x(N,dx)*tau_end+d_1x(N,dx)*tau_end);

vqy_low = zeros(N,1);

heat_low = (1-r).*w;

rad_low = mean(r.*w)*ones(N,1);

%qt_low = qt_low-rad_low;

% top layer: sigma*q = -Uqx -(1-r)w - mean(r(w)w)

q_up = d_2x(N,dx)*phi_end+d_2x(N,dx)*tau_end-tau_end;

qt_up = sigma*q_up;

uqx_up = -(d_3x(N,dx)*phi_end+d_3x(N,dx)*tau_end-d_1x(N,dx)*tau_end);

vqy_up = zeros(N,1);

heat_up = -(1-r).*w;

rad_up = -mean(r.*w)*ones(N,1);

% Residuals

PV_equation_low = (qt_low - uqx_low - vqy_low - heat_low - rad_low);

PV_equation_up = (qt_up - uqx_up - vqy_up - heat_up - rad_up);

% Rescale the different terms and take into account the periodicity

qt_low = [qt_low;qt_low(1)]/max(abs(uqx_low));
vqy_low = [vqy_low;vqy_low(1)]/max(abs(uqx_low));
heat_low = [heat_low;heat_low(1)]/max(abs(uqx_low));
rad_low = [rad_low;rad_low(1)]/max(abs(uqx_low));
Residual_low = [PV_equation_low;PV_equation_low(1)]/max(abs(uqx_low));
uqx_low = [uqx_low;uqx_low(1)]/max(abs(uqx_low));

qt_up = [qt_up;qt_up(1)]/(max(abs(uqx_up)));
vqy_up = [vqy_up;vqy_up(1)]/(max(abs(uqx_up)));
heat_up = [heat_up;heat_up(1)]/(max(abs(uqx_up)));
rad_up = [rad_up;rad_up(1)]/max(abs(uqx_up));
Residual_up = [PV_equation_up;PV_equation_up(1)]/(max(abs(uqx_up)));
uqx_up = [uqx_up;uqx_up(1)]/(max(abs(uqx_up)));

% meridional winds: total and contributions from q1 and q2

% total

v1 = d_1x(N,dx)*(phi_end+tau_end);
v2 = d_1x(N,dx)*(phi_end-tau_end);

% contributions
q1 = d_2x(N,dx)*(phi_end+tau_end)-tau_end;
q2 = d_2x(N,dx)*(phi_end-tau_end)+tau_end;
A2 = d_2x(N,dx); A2(1,:) = 1; A2inv = inv(A2);

phi1 = A2inv*(q1/2);
tau1 = Omega_Solver(q1/2,1,N,dx);

v11 = d_1x(N,dx)*(phi1+tau1);
v12 = d_1x(N,dx)*(phi1-tau1);

phi2 = A2inv*(q2/2);
tau2 = Omega_Solver(-q2/2,1,N,dx);

v21 = d_1x(N,dx)*(phi2+tau2);
v22 = d_1x(N,dx)*(phi2-tau2);


v11 = [v11;v11(1)]/max(abs(v1));
v21 = [v21;v21(1)]/max(abs(v1));

v2 = [v2;v2(1)]/max(abs(v1));
v22 = [v22;v22(1)]/max(abs(v1));
v12 = [v12;v12(1)]/max(abs(v1));
v1 = [v1;v1(1)]/max(abs(v1));

% Plot w and PV-Tendencies of the DRV example 

fig = figure(1);
%x0=10; y0=10; width=1200; height=300; 

x0=10; y0=10; width= 8/7*1300; height= 8/7*300;
set(gcf,'position',[x0,y0,width,height])

% left_color = [0 0 1];
% right_color = [1 0 0];
% set(fig,'defaultAxesColorOrder',[left_color; right_color]);


w = [w;w(1)]/norm(w);

%s1 = subplot(1,3,1);
margx = 0.10;
margy = 0.15;
s1 = subplot_tight(1,3,1,[margy margx]);
%s1 = subplot(1,3,1);

plot(x,w,'linewidth',1.6,'Color','k'); hold on;
xlabel('x'); ylabel('w')
% ax = gca;
% ax.XTick = [0,2*pi,4*pi,6*pi,8*pi];
% ax.XTickLabel = {'0','2\pi','4\pi','6\pi','8\pi'};
title('\rm Vertical Velocity')
set(gca,'fontsize', 15);
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'linewidth',1.7)
xlim([0 L])
annotation('textbox', [0.058, 0.91, 0, 0], 'string', '(a)','fontsize',15)
ylim([-0.05 0.4])
yticks([0 0.1 0.2 0.3 0.4])
yticklabels({'0','0.1','0.2','0.3','0.4'})
%saveas(gcf,'DRV_Paper/Figures/w_DRV_example','epsc')
%set(s1,'position',get(s1,'position')-[0.03 0 0 0])

% Plot PV in the two Layers Plus Diabatic PV Generation

[val,index] = max(abs(uqx_low));
indR = index+100;
indL = index-100;


%fig = subplot(1,3,2);
subplot_tight(2,3,2,[margy margx]);
%fig = subplot(2,3,2);
%set(fig,'position',get(fig,'position')+[-0.02 0 0 0])


%ax1 = gca;
%ax(1).YAxis(1).LineWidth = 1.5;
%yyaxis(ax1,'right')

plot(x(indL:indR),qt_up(indL:indR),'b','linewidth',1.6); hold on
plot(x(indL:indR),heat_up(indL:indR)+rad_up(indL+indR),'b--','linewidth',1.6); hold on;

xlim([x(indL) x(indR)])
set(gca,'XColor','none')
%ylim([-1.2 mean(heat_low)]);
ylim([-1.2 0.014])
yticks([-1 -0.5 0])
yticklabels({'-1','-0.5','0'})
%ylim([-mean(heat_low) disp+mean(heat_low)]);

list = {'$$q_1$$','$$\dot{q}_{1,\mathrm{diab}}$$'};
%columnlegend(1,cellstr(list),'Position',[0 0.43 0.05 0.2],'interpreter','latex','FontSize',9);
legend(list,'interpreter','latex','location','east','FontSize',15);
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'bottom')

set(gca,'XTick',[]);
ax = gca;
ax.YAxis.LineWidth = 1.7;
set(gca,'fontsize', 15);
%set(gca,'linewidth',1.5)



% H=gca;
% H.LineWidth=1.3; %change to the desired value  
t = title('\rm PV Anomalies & Diabatic Generation');
%set(t,'position',get(t,'position')+[-0.1 0 0])


annotation('textbox', [0.365, 0.91, 0, 0], 'string', '(b)','fontsize',15)

subplot_tight(2,3,5,[margy margx]);
%fig = subplot(2,3,5);
%set(fig,'position',get(fig,'position')+[-0.02 0 0 0])


plot(x(indL:indR),qt_low(indL:indR),'r','linewidth',1.6); hold on;
plot(x(indL:indR),heat_low(indL:indR)+rad_low(indL:indR),'r--','linewidth',1.6); hold on;

xlim([x(indL) x(indR)])
%ylim([-mean(heat_low) 1.2])
ylim([-0.014 1.2]);
yticks([0 0.5 1])
yticklabels({'0','0.5','1'})


list = {'$$q_2$$','$$\dot{q}_{2,\mathrm{diab}}$$'};
%columnlegend(1,cellstr(list),'Position',[0 0.43 0.05 0.2],'interpreter','latex','FontSize',9);
legend(list,'interpreter','latex','location','east','FontSize',15);
legend boxoff
xlabel('x')
ax = gca;
ax.YAxis.LineWidth = 1.7;
set(gca,'FontSize',15)
%set(gca,YAxis,'linewidth',1.5);
%set(gca,'fontsize', 12);
%set(gca,'linewidth',1.5)

set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
%set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')



% plot the meridional winds and their decomposition
% into contributions from q1 and q2

%s3 = subplot(1,3,3);

s3 = subplot_tight(2,3,3,[margy margx]);
%s3 = subplot(2,3,3);

%indL = 200;
%indR = 399;
% indL = 1;
% indR = 1001;



plot(x(indL:indR),v1(indL:indR),'b','linewidth',1.6); hold on;
plot(x(indL:indR),v11(indL:indR),'b--','linewidth',1.6); hold on;
%plot(x(indL:indR),v21(indL:indR),'b:','linewidth',1.6); hold on;
p1 = plot(x(indL:indR),zeros(length(indL:indR),1),'k--','linewidth',1.6); hold on;


set(get(get(p1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%y1 = get(gca,'ylim'); hold on;

set(get(get(p1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

xlim([x(indL) x(indR)])
set(gca,'XTick',[]);
set(gca,'XColor','none')

%list = {'$$v_1$$','$$q_1\rightarrow v_1$$','$$q_2\rightarrow v_1$$','$$v_2$$','$$q_2\rightarrow v_2$$','$$q_1\rightarrow v_2$$'};
list = {'$$v_1$$','$$q_1\rightarrow v_1$$'};
legend(list,'interpreter','latex','Position',[0.86 0.79 0 0],'FontSize',15);
legend boxoff
set(gca,'FontSize',15)
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'linewidth',1.7)
title('\rm Meridional Winds')
yticks([-0.5 0 0.5 1])
yticklabels({'-0.5','0','0.5','1'})

annotation('textbox', [0.665, 0.91, 0, 0], 'string', '(c)','fontsize',15)




s3 = subplot_tight(2,3,6,[margy margx]);
%s3 = subplot(2,3,6);

plot(x(indL:indR),v2(indL:indR),'r','linewidth',1.6); hold on;
plot(x(indL:indR),v22(indL:indR),'r--','linewidth',1.6); hold on;
%plot(x(indL:indR),v12(indL:indR),'r:','linewidth',1.6); hold on;
p1 = plot(x(indL:indR),zeros(length(indL:indR),1),'k--','linewidth',1.6);

[pholder,qmax] = max(heat_low);
%p2 = plot([x(qmax) x(qmax)],y1,'k:','linewidth',1.6);


set(get(get(p1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%set(get(get(p2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');


xlabel('x')
xlim([x(indL) x(indR)])



%list = {'$$v_1$$','$$q_1\rightarrow v_1$$','$$q_2\rightarrow v_1$$','$$v_2$$','$$q_2\rightarrow v_2$$','$$q_1\rightarrow v_2$$'};
list = {'$$v_2$$','$$q_2\rightarrow v_2$$'};
legend(list,'interpreter','latex','Position',[0.86 0.22 0 0],'FontSize',15);
legend boxoff

set(gca,'FontSize',15)
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'linewidth',1.7)
yticks([-0.5 0 0.5 1])
yticklabels({'-0.5','0','0.5','1'})
ylim([-1 1])

% note: on matlab 2017b the linewidth of the y-axis on the second panel
% is not properly rendered on the .eps or .pdf save. Seems to be a bug
% that does not appear on newer matlab versions

saveas(gcf,'DRV_Paper/Figures/DRV_mode_example_three_panel_split_up_down','epsc')

% Plot PV-Budget in the Lower Layer
figure(2)
[val,index] = max(abs(uqx_low));
indR = index+60;
indL = index-60;

%subplot(2,2,3)
plot(x(indL:indR),qt_low(indL:indR),'b','linewidth',1.6); hold on
plot(x(indL:indR),uqx_low(indL:indR),'r--','linewidth',1.6); hold on;
plot(x(indL:indR),heat_low(indL:indR),'k-.','linewidth',1.6); hold on;
%plot(x(indL:indR),rad_low(indL:indR),'k:','linewidth',1.6); hold on
xlabel('x'); ylabel('PV tendency')
xlim([x(indL) x(indR)])
ylim([-1.2 1.2])
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'linewidth',1.5)
%h = legend({'\sigma q_2','q_{2x}','(1-r)w'},'Location','northeast'); %,'Orientation','horizontal');
h = legend({'$$\sigma q_2$$','$$q_{2x}$$','$$(1-r)w$$'},'Location','northeast','interpreter','latex'); %,'Orientation','horizontal');
%h = legend({'$$\sigma q_2$$','$$q_{2x}$$','$$(1-r)w$$','$$\overline{r(w)w}$$'},'Location','northeast'); %,'Orientation','horizontal');
set(gca,'fontsize', 12);
legend boxoff
set(h,'Interpreter','latex')
h.FontSize = 14;
%title('PV tendencies in lower layer')
%annotation('textbox', [0.05, 0.5, 0, 0], 'string', '\bf{(c)}','fontsize',12)
%saveas(gcf,'DRV_Paper/Figures/DRV_mode_example_lower_tendencies','epsc')


% Upper PV-Budget for completion
figure(3)
%subplot(2,2,4)

plot(x(indL:indR),qt_up(indL:indR),'linewidth',1.6); hold on
plot(x(indL:indR),uqx_up(indL:indR),'--','linewidth',1.6); hold on;
plot(x(indL:indR),heat_up(indL:indR),'-.','linewidth',1.6); hold on;
%plot(x(indL:indR),rad_up(indL:indR),'k:','linewidth',1.6); hold on
xlabel('x'); ylabel('PV tendency')
xlim([x(indL) x(indR)])
ylim([-1.2 1.2])
set(gca,'fontsize', 12);
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'linewidth',1.5)
h = legend({'$$\sigma q_1$$','$$-q_{1x}$$','$$-(1-r)w$$'},'Location','northeast');
%h = legend({'$$\sigma q_1$$','$$-q_{1x}$$','$$(1-r)w$$','$$\overline{r(w)w}$$'},'Location','northeast'); %,'Orientation','horizontal');
set(h,'Interpreter','latex')
legend boxoff
title('PV tendencies in upper layer')
%annotation('textbox', [0.5, 0.50, 0, 0], 'string', '\bf{(d)}','fontsize',12)

%saveas(gcf,'DRV_Paper/Figures/DRV_mode_example','epsc')

% % Plot PV in the two Layers Plus Diabatic PV Generation
% 
% [val,index] = max(abs(uqx_low));
% indR = index+100;
% indL = index-100;
% 
% 
% fig = subplot(2,2,2);
% % left_color = [0 0 0];
% % right_color = [0 0 0];
% left_color = [0 0 1];
% right_color = [1 0 0];
% set(fig,'defaultAxesColorOrder',[left_color; right_color]);
% 
% ax1 = gca;
% ax(1).YAxis(1).LineWidth = 1.5;
% yyaxis(ax1,'right')
% 
% plot(x(indL:indR),qt_low(indL:indR),'linewidth',2.5); hold on;
% plot(x(indL:indR),heat_low(indL:indR),'linewidth',2); hold on;
% ylim([min(qt_low(indL:indR)) ...
% max(qt_up(indL:indR)+10*max(abs(qt_low)))])
% 
% ax2 = gca;
% yyaxis(ax2,'left')
% ax(2).YAxis(2).LineWidth = 1.5;
% %plot(x(indL:indR),w(indL:indR)+1.5*max(abs(partial_t_term_low)),'Color','k','linewidth',2.5); hold on
% %plot(x(indL:indR),qt_up(indL:indR)+3*max(abs(qt_low)),'Color','b','linewidth',2.5); hold on
% %plot(x(indL:indR),heat_up(indL:indR)/3+3*max(abs(qt_low)),'Color','b','linewidth',2,'Linestyle',':');
% plot(x(indL:indR),qt_up(indL:indR),'linewidth',2.5); hold on
% plot(x(indL:indR),heat_up(indL:indR),'linewidth',2);
% ylim([-max(qt_up(indL:indR)+10*max(abs(qt_low))) -min(qt_low(indL:indR)) ])
% 
% legend('PV1','Diab. Gen.1','PV2','Diab. Gen.2','Location','East'); legend boxoff
% xlabel('x')
% 
% set(gca,'fontsize', 12);
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% xlim([x(indL) x(indR)])
% % ylim([min(qt_low(indL:indR)) ...
% %     max(qt_up(indL:indR)+3*max(abs(qt_low)))])
% % set(gca,'ytick',[])
% % set(gca,'yticklabel',[])
% title('PV Anomalies + Diabatic Generation')
% %saveas(gcf,'/net/halo/disk28/disk7/mkohl/Mode_2D_3D/AGU21_figures/DRV_mode_2layer','epsc')

figure(100)


plot(x,w,'linewidth',1.6); hold on;
xlabel('x'); ylabel('w')
% ax = gca;
% ax.XTick = [0,2*pi,4*pi,6*pi,8*pi];
% ax.XTickLabel = {'0','2\pi','4\pi','6\pi','8\pi'};
title('Vertical velocity')
set(gca,'fontsize', 12);
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'linewidth',1.5)
xlim([0 L])
%saveas(gcf,'/net/halo/disk28/disk7/mkohl/Mode_2D_3D/AGU21_figures/DRV_mode_2layer_vertical_velocity','epsc')


% yyaxis right
% 
% plot(x(indL:indR),qt_low(indL:indR),'Color','r','linewidth',2.5); hold on;
% plot(x(indL:indR),heat_low(indL:indR),'Color','r','linewidth',2,'linestyle',':'); hold on;
% ylim([min(qt_low(indL:indR)) ...
% max(qt_up(indL:indR)+10*max(abs(qt_low)))])
% yyaxis left
% %plot(x(indL:indR),w(indL:indR)+1.5*max(abs(partial_t_term_low)),'Color','k','linewidth',2.5); hold on
% %plot(x(indL:indR),qt_up(indL:indR)+3*max(abs(qt_low)),'Color','b','linewidth',2.5); hold on
% %plot(x(indL:indR),heat_up(indL:indR)/3+3*max(abs(qt_low)),'Color','b','linewidth',2,'Linestyle',':');
% plot(x(indL:indR),qt_up(indL:indR),'Color','b','linewidth',2.5); hold on
% plot(x(indL:indR),heat_up(indL:indR),'Color','b','linewidth',2,'Linestyle',':');
% ylim([-max(qt_up(indL:indR)+10*max(abs(qt_low))) -min(qt_low(indL:indR)) ])
% 
% legend('PV1','Diab. Gen.1','PV2','Diab. Gen.2','Location','East'); legend boxoff
% xlabel('x')
% set(gca,'fontsize', 12);
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% xlim([x(indL) x(indR)])
% % ylim([min(qt_low(indL:indR)) ...
% %     max(qt_up(indL:indR)+3*max(abs(qt_low)))])
% % set(gca,'ytick',[])
% % set(gca,'yticklabel',[])
% title('DRV 2-Layer Mode')
%saveas(gcf,'DRV_Paper/Figures/DRV_mode_example_pv','epsc')

% plot the meridional winds and their decomposition
% into contributions from q1 and q2

% total

v1 = d_1x(N,dx)*(phi_end+tau_end);
v2 = d_1x(N,dx)*(phi_end-tau_end);

% contributions
q1 = d_2x(N,dx)*(phi_end+tau_end)-tau_end;
q2 = d_2x(N,dx)*(phi_end-tau_end)+tau_end;
A2 = d_2x(N,dx); A2(1,:) = 1; A2inv = inv(A2);

phi1 = A2inv*(q1/2);
tau1 = Omega_Solver(q1/2,1,N,dx);

v11 = d_1x(N,dx)*(phi1+tau1);
v12 = d_1x(N,dx)*(phi1-tau1);

phi2 = A2inv*(q2/2);
tau2 = Omega_Solver(-q2/2,1,N,dx);

v21 = d_1x(N,dx)*(phi2+tau2);
v22 = d_1x(N,dx)*(phi2-tau2);


v11 = [v11;v11(1)]/max(abs(v1));
v21 = [v21;v21(1)]/max(abs(v1));

v2 = [v2;v2(1)]/max(abs(v1));
v22 = [v22;v22(1)]/max(abs(v1));
v12 = [v12;v12(1)]/max(abs(v1));
v1 = [v1;v1(1)]/max(abs(v1));

figure(200)

disp = 2;
%indL = 200;
%indR = 399;
indL = 1;
indR = 1001;

% ax2 = gca;
% yyaxis(ax2,'left')
% ax(2).YAxis(2).LineWidth = 1.5;
plot(x(indL:indR),v1(indL:indR)+disp,'b','linewidth',1.6); hold on;
plot(x(indL:indR),v11(indL:indR)+disp,'b--','linewidth',1.6); hold on;
plot(x(indL:indR),v21(indL:indR)+disp,'b:','linewidth',1.6); hold on;
%plot(x,[v11;v11(1)]+[v21;v21(1)],'k--');
% 
% ax1 = gca;
% ax(1).YAxis(1).LineWidth = 1.5;
% yyaxis(ax1,'right')
plot(x(indL:indR),v2(indL:indR),'r','linewidth',1.6); hold on;
plot(x(indL:indR),v22(indL:indR),'r--','linewidth',1.6); hold on;
plot(x(indL:indR),v12(indL:indR),'r:','linewidth',1.6); hold on;
%plot(x,[v22;v22(1)]+[v12;v12(1)],'k')
xlabel('x')
xlim([x(indL) x(indR)])
legend('v1','q1->v1','q2->v1','v2','q2->v2','q1->v2');
legend('Location','East')
legend boxoff
set(gca,'FontSize',12)
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'linewidth',1.5)
title('\rm Meridional Winds')
yticks([-0.5 0 0.5 1 1.5 2 2.5 3])
yticklabels({'-0.5','0','0.5','1','-0.5','0','0.5','1'})

% Plot w for the MM

load('DRV_Paper/MM_example.mat')

w = [w;w(1)]/norm(w);

figure(5)
plot(x,w,'linewidth',1.6); hold on;
xlabel('x'); ylabel('w')
% ax = gca;
% ax.XTick = [0,2*pi,4*pi,6*pi,8*pi];
% ax.XTickLabel = {'0','2\pi','4\pi','6\pi','8\pi'};
title('Untilted Model')
set(gca,'fontsize', 12);
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'linewidth',1.5)
xlim([0 L])
%saveas(gcf,'DRV_Paper/Figures/w_MM_example','epsc')



% make a plot of diabatic heating from one vs. two anomalies

% RHS = 2*d_3x(N,dx)*phi2-d_1x(N,dx)*phi2;
% 
% w2 = Omega_Solver(RHS,0.01,N,dx);
% 
% r2 = r_factor(w2,0.01);
% 
% heat_low2 = (1-r2).*w2;
% 
% heat_low2 = [heat_low2; heat_low2(1)]/(max(abs(heat_low2)));
% 
% figure(1000)
% plot(x(indL:indR),qt_low(indL:indR),'r','linewidth',1.6); hold on;
% plot(x(indL:indR),heat_low(indL:indR),'b--','linewidth',1.6); hold on;
% plot(x(indL:indR),heat_low2(indL:indR),'k--','linewidth',1.6);
% xlim([x(indL) x(indR)])
% xlabel('x'); ylabel('w');
% legend('q2','full diab','q_2 diab')
% legend boxoff


