% Solve for the coefficients in the DRV-Model
close all;
clear;
%v = VideoWriter('/disk7/mkohl/Mode_Zurita_Gotor/Plots_Proposal_Main/Large_Domain/DRV_model.avi','Motion JPEG AVI');
%v.Quality = 100;
%v.FrameRate = 1;
%open(v)

% Plot Growthrate and Asymmetry vs. r

% define r-factor

R = [1.0,0.8,0.6,0.4,0.3,0.2,0.1,0.05,0.02,0.01];

% define q

q = 1;

% define variables

growth_rate = zeros(length(R),1); lambda = zeros(length(R),1); 

for ii = 10:10

%x0 = [0.1,0.2,1,1];

x0 = [0.1,0,1,1];

x = fsolve(@(x) root_solver(x,R(ii),q),x0);
s = x(1); b = x(2); c1 = x(3); c2 = x(4);

k1 = 1/sqrt(2*R(ii))*sqrt(R(ii)*s^2+sqrt((R(ii)*(s^2+2)-1)^2-4*R(ii)*s^2)+2*R(ii)-1);

k2 = 1/sqrt(2*R(ii))*sqrt(R(ii)*s^2-sqrt((R(ii)*(s^2+2)-1)^2-4*R(ii)*s^2)+2*R(ii)-1);

x1 = linspace(-10,0,100);
x2 = linspace(0,x(2),100);
x3 = linspace(x(2),10,100);

w1 = q*s/(s^2-1)*(exp(s*x1)-exp(x1));
w2 = c1*cosh(k1*(x2-b/2))+c2*cosh(k2*(x2-b/2));
w3 = -q*s/(s^2-1)*(exp(-(x3-b))-exp(-s*(x3-b)));

lambda1 = -exp(-s*b)/(k1^2-s^2)*(c1*k1*sinh(k1*b/2)*(R(ii)+1)+c1*s*cosh(k1*b/2)*(R(ii)-1))...
-exp(-s*b)/(k2^2-s^2)*(c2*k2*sinh(k2*b/2)*(R(ii)-1)+c2*s*cosh(k2*b/2)*(R(ii)-1));

q2_M = lambda1*exp(s*x2)+1/(k1^2-s^2)*(c1*k1*sinh(k1*(x2-b/2))*(R(ii)+1)+c1*s*cosh(k1*(x2-b/2))*(R(ii)-1))...
+1/(k2^2-s^2)*(c2*k2*sinh(k2*(x2-b/2))*(R(ii)-1)+c2*s*cosh(k2*(x2-b/2))*(R(ii)-1));

q2_L = q*exp(s*x1);

q2_R = zeros(size(x3));

q2_theory = [q2_L,q2_M,q2_R];

x_axis = [x1,x2,x3];
w_theory = [w1,w2,w3];

figure(1)
plot(x_axis,real(w_theory)); 
xlabel('x')
ylabel('w')
title(['r=',num2str(R(ii))])
pause(2)
%frame = getframe(gcf);
%writeVideo(v,frame);

growth_rate(ii) = x(1);
lambda(ii) = asymmetry(real(w_theory));

end
%close(v);

% Plot Growth Rate

s_DRV = [0,0.3439,0.5371, 0.7017, 0.7860, 0.8836, 1.0222, 1.1395, ...
    1.2690,1.3510];

figure(1)
semilogx(R,s_DRV,'Color','blue','linewidth',1.6); hold on;
semilogx(R,growth_rate,'Color','red','linewidth',1.6); 
xlabel('Reduction factor r','interpreter','latex')
ylabel('\sigma','interpreter','tex')
title('Growthrate')
legend('Numerical','Theory')
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);
set(gca,'linewidth',1.5)


% Plot Asymmetry

lambda_DRV = [0.5,0.7442,0.8272,0.8895,0.9199,0.9478,0.9687,0.9766,...
 0.9829,0.9885];

figure(2)
semilogx(R,lambda_DRV,'Color','blue','linewidth',1.6); hold on;
semilogx(R,lambda,'Color','red','linewidth',1.6); 
xlabel('Reduction factor r','interpreter','latex')
ylabel('Asymmetry parameter $\lambda$','interpreter','latex')
title('Asymmetry')
legend('Numerical','Theory')
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);
set(gca,'linewidth',1.5)
ylim([0.5,1]) 


% Plot Numerical vs. Theoretical Profile for r = 0.01

% Define Grid for DRV-vortex Simulation

load('DRV.mat','w')

N = 1000; 
L = 8*pi; 
dx = L/N;
w = [w;w(1)];
indexL = 200; indexR = 200;

% Define Constants
q = 1;
r = 0.01;

% Nonlinear Solution

x0 = [2,0.2,1,1];
x = fsolve(@(x) root_solver(x,r,q),x0);
sigma = x(1); b = x(2); c1 = x(3); c2 = x(4);

k1 = 1/sqrt(2*r)*sqrt(r*sigma^2+sqrt((r*(sigma^2+2)-1)^2-4*r*sigma^2)+2*r-1);

k2 = 1/sqrt(2*r)*sqrt(r*sigma^2-sqrt((r*(sigma^2+2)-1)^2-4*r*sigma^2)+2*r-1);

x1 = linspace(-indexL*dx,0,indexL+1);
x2 = linspace(0,x(2),16);
x3 = linspace(x(2),x(2)+indexR*dx,indexR+1);
xx = [x1,x2,x3];

w1 = q*sigma/(sigma^2-1)*(exp(sigma*x1)-exp(x1));

w2 = c1*cosh(k1*(x2-b/2))+c2*cosh(k2*(x2-b/2));

w3 = -q*sigma/(sigma^2-1)*(exp(-(x3-b))-exp(-sigma*(x3-b)));

w_theory = [w1,w2,w3];

figure(3)
plot(xx,w(271-(indexL+1):286+(indexR+1))/norm(w(271-(indexL+1):286+(indexR+1))),'Color','red','linewidth',1.6); hold on;
plot(xx,w_theory/norm(w_theory),'--','Color','blue','linewidth',1.6); 
xlim([xx(1) xx(end)])
xlabel('x')
ylabel('w')
legend('Numerical','Theory')
set(gca,'fontsize', 12);
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
legend boxoff
%saveas(gcf,'/disk7/mkohl/Mode_Zurita_Gotor/Plots_Proposal_Main/Large_Domain/DRV_model','epsc')









