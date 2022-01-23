% Compare the Time-Marching Approach to the Theoretical Prediction for the
% DRV-Example

clear; close all;

% Load in DRV-example

load('/net/halo/disk28/disk7/mkohl/Mode_Kohl_OGorman/DRV_Paper/DRV_example.mat')

w = w/norm(w); r = R;

%% Calculate the Theoretical Prediction

% Use asymptotic solution as the initial guess
c0 = 1/2*(1+sqrt(5)); c1 = pi/2*(1+c0)*(c0^2-1)/(1-2*c0);
d2 = c1+pi*(1+2*c0^2)/4;

s_guess = c0+c1*sqrt(r);
b_guess = pi/2*sqrt(r)+(1+c0)*r+d2*r.^(3/2);

% Solve for sigma and b using the two-tangent equations:
x0 = [s_guess,b_guess];
sol = fsolve(@(x) root_solver_3(x,r),x0);

s_theory = sol(1); b_theory = sol(2); 

% Calculate the Wavenumbers k1,k2 and Amplitudes c1 and c2

[k1,k2] = Wavenumbers(s_theory,r);
c1 = 1/(cos(k1*b_theory)*(k2^2-k1^2))*(s_theory^2-1)/R;
c2 = 1/(cos(k2*b_theory)*(k2^2-k1^2))*(1-s_theory^2)/R;

%% Make the Plot for Comparisons

% Define the x-axis

index = find(w>0); 

xu = linspace(-b_theory,b_theory,length(index)+2); 
xd = linspace(b_theory+dx,b_theory+dx+250*dx,252); 

xx = [-flip(xd),xu,xd];

% Define the w-profile based on theoretical solutions

wu = c1*cos(k1*xu)+c2*cos(k2*xu);
wd = (exp(-s_theory*(xd-b_theory))-exp(-(xd-b_theory)));

ww = [flip(wd),wu,wd];

% Make Plot

figure(1)
plot(xx,w(index(1)-1-252:index(end)+1+252),'Color','red','linewidth',1.6); hold on;
plot(xx,ww/norm(ww),'--','Color','blue','linewidth',1.6); 
xlabel('x');
ylabel('w');
legend('Numerical','Theory')
set(gca,'fontsize', 12);
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
legend boxoff
xlim([-max(xd),max(xd)])
