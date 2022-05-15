% Plot the Analytical Solutions and Dispersion Curve for Paper
close all; clear; 

%load('DRV_Analytical/Time_Marching_L32pi_N1000.mat')
load('DRV_Analytical/Time_Marching_L32pi_N1000_detail2.mat')
%load('DRV_Analytical/Time_Marching_N1000_L32pi_new.mat')
%load('DRV_Analytical/Time_Marching_N2000_16pi.mat')
L = 32*pi;
%R = R(1:end-1); s_tm = sigma(1:end-1); b_tm = b(1:end-1);
R = R(2:end); s_tm = sigma(2:end); b_tm = b(2:end);
%R = R; s_tm = sigma; b_tm = b;
%R = R(5:end-2); s_tm = sigma(5:end-2); b_tm = b(5:end-2);

% Define Properties of fsolve

options = optimoptions('fsolve'); 
options.MaxIterations = 100000;
options.MaxFunctionEvaluations = 100000;

% Reducation factor
r = R;

%Pre-initialize growth rate and half-ascent area

s_disp = zeros(length(r),1); b_disp = zeros(length(r),1);
C1 = zeros(length(r),1); C2 = zeros(length(r),1); a = zeros(length(r),1);
d = zeros(length(r),1);

% Theory: Parameters
c0 = 1/2*(1+sqrt(5)); c1 = pi/2*(1+c0)*(c0^2-1)/(1-2*c0);
d2 = c1+pi*(1+2*c0^2)/4;
c2 = (pi/2*(c0^2*c1+2*c0*c1)+(1+c0)^2*(c0^2-1)+c1^2)/(1-2*c0);

s_asym = c0+c1*sqrt(r);%+c2*r;
b_asym = pi/2*sqrt(r)+(1+c0)*r+d2*r.^(3/2);



for ii = 1:length(R)
    
% Solve the two tangent equations using the asymptotic solution as initial
% guess
    
x0 = [s_asym(ii)+c2*r(ii),b_asym(ii)];
%x0 = [s_tm(ii),b_tm(ii)];
sol = fsolve(@(x) root_solver_3(x,r(ii)),x0,options);
 
s_disp(ii) = sol(1); b_disp(ii) = sol(2); 

% Infer the Wavenumbers and remaining Amplitudes of the Solutions

[k1,k2] = Wavenumbers(s_disp(ii),r(ii));

% Calculate the Amplitudes of the Ascending solution, take into account
% that d needs to be negative when sigma<1 in order for w_down to be
% negative

if s_disp(ii)>=1
    d(ii) = 1;
else
    d(ii) = -1;
end

C1(ii) = d(ii)/(cos(k1*b_disp(ii))*(k2^2-k1^2))*(s_disp(ii)^2-1)/r(ii);
C2(ii) = d(ii)/(cos(k2*b_disp(ii))*(k2^2-k1^2))*(1-s_disp(ii)^2)/r(ii);

% Calculate the a factor

a(ii) = 1/L*(s_disp(ii)*(s_disp(ii)-1)...
    -s_disp(ii)^3*(s_disp(ii)-1)/(s_disp(ii)^2+r(ii)-1));

end

figure(1)
loglog(r,s_asym,'Linewidth',1.5); hold on;
loglog(r,s_disp,'Linewidth',1.5); hold on;
loglog(R,s_tm,'Linewidth',1.5); hold on;
xlabel('r')
ylabel('\sigma')
xlim([min(r) max(r)])
legend('Asymptotic','Dispersion','Time Marching','Location','Northwest')
legend boxoff
set(gca,'fontsize', 12);
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
%saveas(gcf,'DRV_Paper/Figures/sigma_dispersion','epsc')

figure(2)
semilogx(r,b_asym,'Linewidth',1.5); hold on;
semilogx(r,b_disp,'Linewidth',1.5); hold on;
semilogx(R,b_tm,'Linewidth',1.5);
legend('Asymptotic','Dispersion','Time Marching')
legend boxoff
xlabel('r')
ylabel('b')
xlim([min(r) max(r)])
set(gca,'fontsize', 12);
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
%saveas(gcf,'DRV_Paper/Figures/b_dispersion','epsc')

for ii = 1:length(R)

% Define x-axis

%xu = linspace(0,b_disp(ii),100); xd = linspace(b_disp(ii),L,100);

dxx = b_disp(ii)/1000;

xu = 0:dxx:b_disp(ii); xd = b_disp(ii)+dxx:dxx:L;

xx = [-flip(xd),-flip(xu),xu,xd];

% Define w-solutions
[k1,k2] = Wavenumbers(s_disp(ii),R(ii));

wu = a(ii)/(s_disp(ii)^2+R(ii)-1)+C1(ii)*cos(k1*xu)+C2(ii)*cos(k2*xu);
wd = a(ii)/s_disp(ii)^2*(1-exp(-(xd-b_disp(ii))))+d(ii)*(exp(-s_disp(ii)*(xd-b_disp(ii)))-exp(-(xd-b_disp(ii))));

%wu = C1(ii)*cos(k1*xu)+C2(ii)*cos(k2*xu);
%wd = d(ii)*(exp(-s_disp(ii)*(xd-b_disp(ii)))-exp(-(xd-b_disp(ii))));

ww = [flip(wd),flip(wu),wu,wd];

mass(ii) = trapz(xu,wu)+trapz(xd,wd);

% Calculate Mass using matlab's integral function

Wu = @(xx) a(ii)/(s_disp(ii)^2+R(ii)-1)+C1(ii)*cos(k1*xx)+C2(ii)*cos(k2*xx);

Wd = @(xx) a(ii)/s_disp(ii)^2*(1-exp(-(xx-b_disp(ii))))+d(ii)*(exp(-s_disp(ii)*(xx-b_disp(ii)))-exp(-(xx-b_disp(ii))));


mass2(ii) = integral(Wu,0,b_disp(ii))+integral(Wd,b_disp(ii),L);

% 
al(ii) = -(1-1/s_disp(ii)-s_disp(ii)*(s_disp(ii)-1)/(s_disp(ii)^2+r(ii)-1));

% Make plot of the solution

figure(ii+10)
plot(xx,ww,'Linewidth',1.5); 
xlabel('x'); ylabel('w')
%xlim([xx(1) xx(end)])
title(['R=',num2str(R(ii))])

end

% Check the RHS of the tangent equations
[K1,K2,RHS1,RHS2,term1,term2] = deal(zeros(length(R),1));
for ii = 1:length(R)
    
%x(1) = s_tm(ii); x(2) = b_tm(ii); 

%[k1,k2] = Wavenumbers(s_tm(ii),r(ii)); 

x(1) = s_disp(ii); x(2) = b_disp(ii); 

[k1,k2] = Wavenumbers(s_disp(ii),r(ii)); 

K1(ii) = k1; K2(ii) = k2;

Root_1(ii) = 1-R(ii)*(2+s_disp(ii)^2);

Root_2(ii) = sqrt(1+R(ii)^2*(s_disp(ii)^4+4*s_disp(ii)^2)-6*R(ii)*s_disp(ii)^2);

Root_inner(ii) = 1+R(ii)^2*(s_disp(ii)^4+4*s_disp(ii)^2)-6*R(ii)*s_disp(ii)^2;

if isreal(k2)==1

RHS1(ii) =    -r(ii)*k1*k2/(1+x(1))*(1/(r(ii)*k2)-x(1)*k1/(x(1)^2+r(ii)-1));
    

RHS2(ii) = r(ii)*k1*k2/(1+x(1))*(-1/(r(ii)*k1)+x(1)*k1/(x(1)^2+r(ii)-1));

%RHS2(ii) = (-1/(r(ii)*k1)+x(1)*k1/(x(1)^2+r(ii)-1));

term1(ii) = -1/(r(ii)*k1); term2(ii) = x(1)*k1/(x(1)^2+r(ii)-1);

else

[k1,k2] = Wavenumbers_Complex(s_disp(ii),r(ii)); 
    
RHS1(ii) =    -r(ii)*k1*k2/(1+x(1))*(1/(r(ii)*k2)+x(1)*k1/(x(1)^2+r(ii)-1));

%RHS2(ii) = (-1/(r(ii)*k1)+x(1)*k1/(x(1)^2+r(ii)-1));

RHS2(ii) = r(ii)*k1*k2/(1+x(1))*(-1/(r(ii)*k1)+x(1)*k1/(x(1)^2+r(ii)-1));

term1(ii) = -1/(r(ii)*k1); term2(ii) = x(1)*k1/(x(1)^2+r(ii)-1);

end

end

C1_tm = 1./(cos(K1.*b_tm).*(K2.^2-K1.^2)).*(s_tm.^2-1)./r';
C2_tm = 1./(cos(K2.*b_tm).*(K2.^2-K1.^2)).*(1-s_tm.^2)./r';

figure(10)
plot(R,RHS1,'Linewidth',1.5); hold on;
plot(R,ones(size(R)),'r','Linewidth',1.5); hold on;
plot(R,-ones(size(R)),'r','Linewidth',1.5)
xlabel('R')
ylabel('Amp')
title('RHS Tangent Equation2')

figure(100)
plot(R,RHS2,'Linewidth',1.5); hold on;
plot(R,ones(size(R)),'r','Linewidth',1.5); hold on;
plot(R,-ones(size(R)),'r','Linewidth',1.5)
xlabel('R')
ylabel('Amp')
title('RHS Tangent Equation2')


figure(20)
plot(R,real(K2),'Linewidth',1.5); hold on; 
title('K2')
legend('Real'); legend boxoff
xlabel('R')
ylabel('Amp')

figure(21)
plot(R,imag(K2),'Linewidth',1.5); 
title('K2')
legend('Imag'); legend boxoff
xlabel('R')
ylabel('Amp')

figure(22)
plot(R,real(K1),'Linewidth',1.5); hold on; 
title('K1')
legend('Real'); legend boxoff
xlabel('R')
ylabel('Amp')

figure(23)
plot(R,imag(K1),'Linewidth',1.5); 
title('K1')
legend('Imag'); legend boxoff
xlabel('R')
ylabel('Amp')

figure(24)
plot(R,real(Root_inner),'Linewidth',1.5); 
title('Inner Root')
legend('Real Part')
xlabel('R')
ylabel('Amp')

figure(25)
plot(R,imag(Root_inner),'Linewidth',1.5); 
legend('Imag Part')
title('Inner Root')
xlabel('R')
ylabel('Amp')

figure(26)
plot(R,Root_inner); hold on;
plot(R,Root_1)

