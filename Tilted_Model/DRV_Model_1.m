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

% test

x0 = [1.5,0.1];

%x0 = [0.9825,0.9];

%x = fsolve(@(x) root_solver_2(x,0.01),x0);
x(1) = 1.3189; x(2) =  0.1885;
r = 0.01;

k1 = 1/sqrt(2*r)*sqrt(1-r*(x(1)^2+2)+sqrt(1+r^2*(x(1)^4+4*x(1)^2)-6*r*x(1)^2));

k2 = 1/sqrt(2*r)*sqrt(1-r*(x(1)^2+2)-sqrt(1+r^2*(x(1)^4+4*x(1)^2)-6*r*x(1)^2));

F1 = tan(k1*x(2))+r*k1*k2/(x(1)+1)*(1/(r*k2)-k2/x(1));

F2 = tan(k2*x(2))-r*k1*k2/(x(1)+1)*(-1/(r*k1)+k1/x(1));

for ii = 10:10

%x0 = [0.1,0.2,1,1];

x0 = [1,0,1,1];

x = fsolve(@(x) root_solver_1(x,R(ii)),x0);
s = x(1); b = x(2); c1 = x(3); c2 = x(4);

k1 = 1/sqrt(2*R(ii))*sqrt(1-R(ii)*(s^2+2)+sqrt(1+R(ii)^2*(s^4+4*s^2)-6*R(ii)*s^2));

k2 = 1/sqrt(2*R(ii))*sqrt(1-R(ii)*(s^2+2)-sqrt(1+R(ii)^2*(s^4+4*s^2)-6*R(ii)*s^2));


x1 = linspace(0,b,100); x2 = linspace(b,10,100);
x_axis = [-flip(x2),-flip(x1),x1,x2];

w1 = c1*cos(k1*x1)+c2*cos(k2*x1);
w2 = exp(-s*(x2-b))-exp(-(x2-b));

w_theory = [flip(w2),flip(w1),w1,w2];
w_theory = w_theory/norm(w_theory);

figure(1)
plot(x_axis,w_theory); 
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



