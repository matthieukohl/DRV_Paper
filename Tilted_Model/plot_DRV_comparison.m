% Compare the Time-Marching Approach to the Theoretical Prediction for the
% DRV-Example

clear; close all;

% Load in DRV-example

load('/net/halo/disk28/disk7/mkohl/Mode_Kohl_OGorman/DRV_Paper/DRV_example.mat')

w = w/max(w); r = R; 

% Specify the domain from time-marching tobe compared to

[val,xmax] = max(w); Nindex = 170;

LL = Nindex*dx; x = -LL:dx:LL; index = xmax-Nindex:1:xmax+Nindex;

%% Calculate the Theoretical Prediction

% Use asymptotic solution as the initial guess
c0 = 1/2*(1+sqrt(5)); c1 = pi/2*(1+c0)*(c0^2-1)/(1-2*c0);
d2 = c1+pi*(1+2*c0^2)/4;

s_guess = c0+c1*sqrt(r);
b_guess = pi/2*sqrt(r)+(1+c0)*r+d2*r.^(3/2);

x0 = [s_guess,b_guess];

% Solve for sigma and b using the two-tangent equations:

% Define Solver Options
options = optimset('TolX',1e-12,'TolFun',1e-12);

% fminsearch method

%fun = @(x)residual_fn(x,r);
%[sol,fval] = fminsearch(fun,x0,options);

% structure of unphysical second branch solution (w<0);

% x0 = [1.2374,0.5066];
% fun = @(x)residual_fn(x,r);
% [sol,fval] = fminsearch(fun,x0,options);

% fsolve method

sol = fsolve(@(x) root_solver_3(x,r),x0,options);


s_theory = sol(1); b_theory = sol(2); 

% Calculate the Wavenumbers k1,k2 and Amplitudes c1 and c2

[k1,k2] = Wavenumbers(s_theory,r);
c1 = 1/(cos(k1*b_theory)*(k2^2-k1^2))*(s_theory^2-1)/R;
c2 = 1/(cos(k2*b_theory)*(k2^2-k1^2))*(1-s_theory^2)/R;

%% Make the Plot for Comparisons

% Define the x-axis
Nu = 100; Nd = 200;
 
xu = linspace(-b_theory,b_theory,Nu); 
xd = linspace(b_theory,LL,Nd); 

xx = [-flip(xd),xu,xd];

% Define the w-profile based on theoretical solutions

wu = c1*cos(k1*xu)+c2*cos(k2*xu);
wd = (exp(-s_theory*(xd-b_theory))-exp(-(xd-b_theory)));

ww = [flip(wd),wu,wd];

% Make Plot

figure(1)
plot(x,w(index),'r','linewidth',1.6); hold on;
plot(xx,ww/max(ww),'b--','linewidth',1.6); hold on;
xlabel('x');
ylabel('w');
legend('Finite domain','Infinite domain')
set(gca,'fontsize', 12);
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
legend boxoff
xlim([-max(xd),max(xd)])
%saveas(gcf,'DRV_Paper/Figures/DRV_comparison','epsc')

% % Compare theory to DRV-GCM

close all; clear;
load('/net/halo/disk28/disk7/mkohl/Mode_Kohl_OGorman/DRV_Analytical/DRV_GCM.mat')
w = squeeze(-mean(omega(latm,:,5:18,1),3)); w= w/max(w); r = rr; R =r; 
[val,lonm] = max(w);
R_Earth = 6371000.0; 
dx = R_Earth*cosd(lat(latm))*abs(lon(2)-lon(1))*pi/180; dx = dx/LD;

LL = 60*dx; x = -LL:dx:LL; index = lonm-60:1:lonm+60;

% Use asymptotic solution as the initial guess
c0 = 1/2*(1+sqrt(5)); c1 = pi/2*(1+c0)*(c0^2-1)/(1-2*c0);
d2 = c1+pi*(1+2*c0^2)/4;

s_guess = c0+c1*sqrt(r);
b_guess = pi/2*sqrt(r)+(1+c0)*r+d2*r.^(3/2);

x0 = [s_guess,b_guess];

% Solve for sigma and b using the two-tangent equations:

% Define Solver Options
options = optimset('TolX',1e-12,'TolFun',1e-12);

%x = fsolve(@(x) root_solver_3(x,r),x0);

fun = @(x)residual_fn(x,r);
[sol,fval] = fminsearch(fun,x0,options);

s_theory = sol(1); b_theory = sol(2); 

% Calculate the Wavenumbers k1,k2 and Amplitudes c1 and c2

[k1,k2] = Wavenumbers(s_theory,r);
c1 = 1/(cos(k1*b_theory)*(k2^2-k1^2))*(s_theory^2-1)/R;
c2 = 1/(cos(k2*b_theory)*(k2^2-k1^2))*(1-s_theory^2)/R;

% Define the x-axis

Nu = 100; Nd = 200;

xu = linspace(-b_theory,b_theory,Nu); 
xd = linspace(b_theory,LL,Nd); 

xx = [-flip(xd),xu,xd];

% Define the w-profile based on theoretical solutions

wu = c1*cos(k1*xu)+c2*cos(k2*xu);
wd = (exp(-s_theory*(xd-b_theory))-exp(-(xd-b_theory)));

ww = [flip(wd),wu,wd];

% Make Plot

figure(2)
plot(x,w(index),'Color','red','linewidth',1.6); hold on;
plot(xx,ww/max(ww),'--','Color','blue','linewidth',1.6); 
xlabel('x');
ylabel('w');
legend('GCM','Root-Solving')
set(gca,'fontsize', 12);
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
legend boxoff
xlim([-max(xd),max(xd)])

