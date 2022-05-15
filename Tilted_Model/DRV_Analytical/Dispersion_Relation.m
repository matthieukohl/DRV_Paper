%% Solve the Dispersion Relationship for the DRV:
% tan(k1*b) = -rk1*k2/(sigma+1)*(1/r*k2-sigma*k2/(sigma^2+r-1))
% tan(k2*b) = r*k1*k2/(sigma+1)*(-1/r*k1+sigma*k1/(sigma^2+r-1))

clear; close all;

% Define Properties of fsolve

options = optimoptions('fsolve'); 
options.MaxIterations = 100000;
options.MaxFunctionEvaluations = 100000;

M = 100;

%load('Time_Marching_HR.mat')
%load('Time_Marching_L32pi_N1000.mat')
load('DRV_Analytical/Time_Marching_L32pi_N1000_detail2.mat')
L = 32*pi;
R = R(1:end-1); sigma = sigma(1:end-1); b = b(1:end-1);

% Reducation factor
%r = linspace(1,0.4,M);
r = R;
%Pre-initialize growth rate and half-ascent area

sigma1 = zeros(length(r),1); b1 = zeros(length(r),1);

% Theory: Parameters
c0 = 1/2*(1+sqrt(5)); c1 = pi/2*(1+c0)*(c0^2-1)/(1-2*c0);
d2 = c1+pi*(1+2*c0^2)/4;

s_theory = c0+c1*sqrt(r);
b_theory = pi/2*sqrt(r)+(1+c0)*r+d2*r.^(3/2);

lambda = c0^2*(2*c0^2-1)/(2*sqrt(c0^2-1));

RHS = c0^2*(1-2*c0^2)-c0*(1+2*c0^2)-c0/(c0^2-1)+c0^3*(2*c0^2-1)/(2*(c0^2-1))+...
    lambda-lambda*c0/(c0^2-1);


c2 = (pi/2*(c0^2*c1+2*c0*c1)+(1+c0)^2*(c0^2-1)+c1^2)/(1-2*c0);

%x0 = [c0+c1*sqrt(r(1))+c2*r(1);pi/2*sqrt(r(1))+(1+c0)*r(1)+d2*r(1).^(3/2)];
for ii = 1:length(R)
x0 = [c0+c1*sqrt(r(ii))+c2*r(ii);pi/2*sqrt(r(ii))+(1+c0)*r(ii)+d2*r(ii).^(3/2)];
%x0 = [c0+c1*sqrt(r(ii));pi/2*sqrt(r(ii))+(1+c0)*r(ii)+d2*r(ii).^(3/2)];
%x0 = [sigma(ii),b(ii)];
% %x0 = [rand,rand];
sol = fsolve(@(x) root_solver_3(x,r(ii)),x0,options);
% 
 sigma1(ii) = sol(1); b1(ii) = sol(2); [k1,k2] = Wavenumbers(sigma1(ii),r(ii));
 C1(ii) = 1/(cos(k1*b1(ii))*(k2^2-k1^2))*(sigma1(ii)^2-1)/r(ii);
C2(ii) = 1/(cos(k2*b1(ii))*(k2^2-k1^2))*(1-sigma1(ii)^2)/r(ii);
a(ii) = 0;

%x0 = [sigma(ii),b(ii),1,1];
%x0 = [c0+c1*sqrt(r(ii))+c2*r(ii);pi/2*sqrt(r(ii))+(1+c0)*r(ii)+d2*r(ii).^(3/2);1;1];
%x0 = [rand,rand,rand,rand];
%sol = fsolve(@(x) root_solver_2(x,r(ii)),x0,options);

%sigma1(ii) = sol(1); b1(ii) = sol(2); C1(ii) = sol(3); C2(ii) = sol(4);
%a(ii) = 0;

%x0 = [sigma(ii),b(ii),1,1,-5];
%x0 = [sigma(ii),b(ii),1,1,-0.1];
%x0 = [c0+c1*sqrt(r(ii))+c2*r(ii);pi/2*sqrt(r(ii))+(1+c0)*r(ii)+d2*r(ii).^(3/2);1;1;-5];
%sol = fsolve(@(xx) root_solver(r(ii),L,xx),x0,options);

%sigma1(ii) = sol(1); b1(ii) = sol(2); C1(ii) = sol(3); C2(ii) = sol(4); a(ii) = sol(5);

% use previous solution as new guess
%x0 = [sol(1),sol(2)];

% find accuracy of solution

x = [sigma1(ii),b1(ii)];

[k1,k2] = Wavenumbers(sigma1(ii),r(ii)); 

F(1) = tan(k1*x(2))+r(ii)*k1*k2/(1+x(1))*(1/(r(ii)*k2)-k2*x(1)/(x(1)^2+r(ii)-1));
% 
F(2) = tan(k2*x(2))-r(ii)*k1*k2/(1+x(1))*(-1/(r(ii)*k1)+x(1)*k1/(x(1)^2+r(ii)-1));

% if isreal(k2)==1
%     
% [k1,k2] = Wavenumbers(sigma1(ii),r(ii)); 
% 
% F(1) = tan(k1*x(2))+r(ii)*k1*k2/(1+x(1))*(1/(r(ii)*k2)-k2*x(1)/(x(1)^2+r(ii)-1));
% % 
% F(2) = tan(k2*x(2))-r(ii)*k1*k2/(1+x(1))*(-1/(r(ii)*k1)+x(1)*k1/(x(1)^2+r(ii)-1));
% 
% else
% 
% [k1,k2] = Wavenumbers_Complex(sigma1(ii),r(ii));   
%     
% F(1) = tan(k1*x(2))+r(ii)*k1*k2/(1+x(1))*(1/(r(ii)*k2)+k2*x(1)/(x(1)^2+r(ii)-1));
% % 
% F(2) = tanh(k2*x(2))-r(ii)*k1*k2/(1+x(1))*(-1/(r(ii)*k1)+x(1)*k1/(x(1)^2+r(ii)-1));
% 
% end
    

accuracy1(ii) = F(1); accuracy2(ii) = F(2);
end

% for ii = 1:M
% x0 = [c0+c1*sqrt(r(ii))+c2*r(ii),pi/2*sqrt(r(ii))+(1+c0)*r(ii)+d2*r(ii).^(3/2),0.3061,0.07002];
% sol = fsolve(@(x) root_solver_2(x,r(ii)),x0);
% sigma(ii) = sol(1); b(ii) = sol(2); C1(ii) = sol(3); C2(ii) = sol(4);
% 
% % use previous solution as new guess
% %x0 = [sol(1),sol(2)];
% end

% Time-Marching Solutions

r_time = [0.1000,0.0800,0.0600,0.0400,0.0100,0.0080,0.0060,0.0040,0.0010];
b_time = [0.9048,0.7791,0.6283,0.4524,0.2011,0.1759,0.1508,0.1257,0.0503];
s_time = [1.0419,1.0811,1.1290,1.1914,1.3620,1.3836,1.4092,1.4414,1.5247];

%load('Time_Marching.mat')
% ii =10;
% [k1,k2] = Wavenumbers(sigma(ii),r(ii));
% 
% x = [sigma(ii),b(ii)];
% 
% F(1) = tan(k1*x(2))+r(ii)*k1*k2/(1+x(1))*(1/(r(ii)*k2)-k2*x(1)/(x(1)^2+r(ii)-1));
% % 
% F(2) = tanh(k2*x(2))-r(ii)*k1*k2/(1+x(1))*(-1/(r(ii)*k1)+x(1)*k1/(x(1)^2+r(ii)-1));

figure(1)
semilogx(r,s_theory); hold on;
semilogx(r,sigma1); hold on;
semilogx(R,sigma); hold on;
xlabel('r')
xlim([min(r) max(r)])
legend('Asymptotic','Dispersion','Time Marching')
legend boxoff
title('Growth Rate')

figure(2)
semilogx(r,b_theory); hold on;
semilogx(r,b1); hold on;
semilogx(R,b);
legend('Asymptotic','Dispersion','Time Marching')
legend boxoff
xlabel('r')
xlim([min(r) max(r)])
title('Half-Ascent Area')
% 
% [k1,k2] = Wavenumbers(sigma1,r');
% 
% figure(3)
% plot(r,k1); hold on;
% legend('k1')
% legend boxoff
% title('Wavenumbers')
% 
% figure(4)
% plot(r,k2); hold on;
% legend('k2')
% legend boxoff
% title('Wavenumbers')
% 
% figure(5)
% plot(r,accuracy1); hold on
% legend('F1')
% legend boxoff
% title('Accuracy')
% 
% figure(6)
% plot(r,accuracy2); hold on
% legend('F2')
% legend boxoff
% title('Accuracy')




% Plot w-solution 



for ii = 1:length(R)

% Define x-axis

%xu = linspace(0,b1(ii),100); xd = linspace(b1(ii),L,100);

dxx = b1(ii)/1000;

xu = 0:dxx:b1(ii); xd = b1(ii)+dxx:dxx:L;

xx = [-flip(xd),-flip(xu),xu,xd];
%xx = [xu,xd];

% Define w-solutions
[k1,k2] = Wavenumbers(sigma1(ii),R(ii));

wu = a(ii)/(sigma1(ii)^2+R(ii)-1)+C1(ii)*cos(k1*xu)+C2(ii)*cos(k2*xu);
wd = a(ii)/sigma1(ii)^2*(1-exp(-(xd-b1(ii))))+(exp(-sigma1(ii)*(xd-b1(ii)))-exp(-(xd-b1(ii))));

ww = [flip(wd),flip(wu),wu,wd];
%ww = [wu,wd];

% mass(ii) = a(ii)*b1(ii)/(r(ii)+sigma1(ii)^2-1)+C1(ii)/k1*sin(k1*b1(ii))+...
%     C2(ii)/k2*sin(k2*b1(ii))+a(ii)/(sigma1(ii)^2)*(L-b1(ii))+a(ii)/sigma1(ii)^2*(exp(-(L-b1(ii)))-1)+...
%     -1/sigma1(ii)*(exp(-sigma1(ii)*(L-b1(ii)))-1)+exp(-(L-b1(ii)))-1;
%mass(ii) = C1(ii)/k1*sin(k1*b1(ii))+C2(ii)/k2*sin(k2*b1(ii))+1/sigma1(ii)-1; %trapz(xx,ww);
%mass(ii) = trapz(xx,ww);

mass(ii) = trapz(xu,wu)+trapz(xd,wd);

Wu = @(xx) a(ii)/(sigma1(ii)^2+R(ii)-1)+C1(ii)*cos(k1*xx)+C2(ii)*cos(k2*xx);

Wd = @(xx) a(ii)/sigma1(ii)^2*(1-exp(-(xx-b1(ii))))+(exp(-sigma1(ii)*(xx-b1(ii)))-exp(-(xx-b1(ii))));

mass2(ii) = integral(Wu,0,b1(ii))+integral(Wd,b1(ii),L);
% Make plot of the solution

figure(ii+10)
plot(xx,ww,'Linewidth',1.5); 
xlabel('x'); ylabel('w')
%xlim([xx(1) xx(end)])
title(['R=',num2str(R(ii))])

end

%close all;

% Check the RHS of the tangent equations

for ii = 1:length(R)
    
x(1) = sigma(ii); x(2) = b(ii); 

[k1,k2] = Wavenumbers(sigma(ii),r(ii)); 

if isreal(k2)==1

RHS(ii) = r(ii)*k1*k2/(1+x(1))*(-1/(r(ii)*k1)+x(1)*k1/(x(1)^2+r(ii)-1));

else

[k1,k2] = Wavenumbers_Complex(sigma(ii),r(ii)); 
    
RHS(ii) = r(ii)*k1*k2/(1+x(1))*(-1/(r(ii)*k1)+x(1)*k1/(x(1)^2+r(ii)-1));

end

end












% for ii = 2:2
%     
% x = [sigma1(ii),b1(ii),C1(ii),C2(ii),a(ii)]; r = R(ii);
% 
% 
% k1 = 1/sqrt(2*r)*sqrt(1-r*(2+x(1)^2)+sqrt(1+r^2*(x(1)^4+4*x(1)^2)-6*r*x(1)^2));
% 
% 
% k2 = 1/sqrt(2*r)*sqrt(1-r*(2+x(1)^2)-sqrt(1+r^2*(x(1)^4+4*x(1)^2)-6*r*x(1)^2));
% 
% 
% F(1) = x(5)/(r+x(1)^2-1)+x(3)*cos(k1*x(2))+x(4)*cos(k2*x(2));
% 
% F(2) = -r*x(3)*k1*sin(k1*x(2))-r*x(4)*k2*sin(k2*x(2))-x(5)/x(1)^2-(1-x(1));
% 
% F(3) = x(5)*x(2)/(r+x(1)^2-1)+x(3)/k1*sin(k1*x(2))+x(4)/k2*sin(k2*x(2))+...
%     x(5)/x(1)^2*(L-x(2))+x(5)/x(1)^2*(exp(-(L-x(2)))-1)+...
%     -1/x(1)*(exp(-x(1)*(L-x(2)))-1)+(exp(-(L-x(2)))-1);
% 
% F(4) = -r*x(3)*k1^2*cos(k1*x(2))-r*x(4)*k2^2*cos(k2*x(2))+x(5)/x(1)^2-(x(1)^2-1);
% 
% F(5) = x(5)*x(2)*r/(r+x(1)^2-1)+r*x(3)/k1*sin(k1*x(2))+r*x(4)/k2*sin(k2*x(2))+...
%    x(5)/x(1)^2*(L-x(2))+x(5)/x(1)^2*(exp(-(L-x(2)))-1)+...
%    -1/x(1)*(exp(-(L-x(2)))-1)+(exp(-(L-x(2)))-1)-x(5)*L;
% 
% accuracy(ii) = max(abs(F));
%    
% end
 
%    
% figure(10)
% plot(R,accuracy)

