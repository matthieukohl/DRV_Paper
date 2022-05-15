% Compare the root solving solutions (sigma,b,w) to the dispersion relation 
% for an infinite domain to the time marching solutions on a finite domain


clear; close all;

%load('DRV_Analytical/Time_Marching_L32pi_N1000_detail2.mat')
%load('DRV_Analytical/Time_Marching_N2000_16pi.mat')
%load('DRV_Analytical/Time_Marching_L16pi_N2000.mat')
%load('DRV_Analytical/Time_Marching_L8pi_N1000.mat')

%load('DRV_Paper/Time_Marching_L32pi_N800.mat')

load('DRV_Paper/Time_Marching_L32pi_N1200_more_values.mat')
%load('DRV_Paper/Time_Marching_L8pi_N200_more_r_values.mat')



R = R(1:end); s_tm = sigma(1:end); b_tm = b(1:end); clear b;
L = 32*pi*10; 

% growth rate at r = 1 is negative, so set to zero and exclude its b value

s_tm(end) = NaN; b_tm(end) = NaN;

% Define r

r = linspace(0.001,0.379,4e3);
%r = linspace(0.01,0.38,1000);
%r = linspace(0.001,0.01,1000);
%r = flip(r);
%r = linspace(0.4,0.8,1000);

%r_start = (3-sqrt(5))/2;
%r = linspace(r_start,0.8,1000);%
%r = flip(r);

% Define Solver Options
%options = optimset('TolX',1e-12,'TolFun',1e-12);
options = optimoptions('fsolve'); 
%options.MaxIterations = 100000;
%options.MaxFunctionEvaluations = 100000;

% Define Initial Guess

%x0(1) = 1.5306;
%x0(2) = 0.0526;
x0(1) = 1.53;
x0(2) = 0.06;

%x0(1) = 1.3; 
%x0(2) = 0.2; % second solution r = 0.001

%x0(1) = 0.3;
%x0(2) = -5;
%x0(1) = -2;
%x0(2) = 20;

%x0(1) = 0;
%x0(2) = -30;

%x0(1) = 0.79;
%x0(2) = 4;

% Pre-define values

[s,b,fval,C1,C2,d,a,mass,mass2] = deal(zeros(1,length(r)));

% asymptotic predictions: sigma = c0, b = pi/2*sqrt(r)+(1+c0)*r

c0 = 1/2*(1+sqrt(5)); b_asym = pi/2*sqrt(r)+(1+c0)*r;

 
 for i = 1:length(r)

% linear interpolation as initial guess


if i>2
x0(1) = 2*s(i-1)-s(i-2); 
x0(2) = 2*b(i-1)-b(i-2); 
end

% Solve the tan-equations for growth rate s and half-ascent area b

if i==400
    test =0;
end

%fun = @(x)residual_fn(x,r(i));
%[x,fval(i)] = fminsearch(fun,x0,options);
sol = fsolve(@(x) root_solver_3(x,r(i)),x0,options);
%sol = fsolve(@(x) root_solver_beyond_cut_off(x,r(i)),x0,options);


s(i) = sol(1);
b(i) = sol(2);


% calculate and plot the w-profiles (optional)

% Calculate various coefficients in the solutions

[k1,k2] = Wavenumbers(s(i),r(i));

% sign switch for d needed to ensure that w<0 in the descending branch 
% when sigma<1

if s(i)>=1
    d(i) = 1;
else
    d(i) = -1;
end

C1(i) = d(i)/(cos(k1*b(i))*(k2^2-k1^2))*(s(i)^2-1)/r(i);
C2(i) = d(i)/(cos(k2*b(i))*(k2^2-k1^2))*(1-s(i)^2)/r(i);
a(i) = d(i)/L*((r(i)-1)*s(i)*(s(i)-1)/(r(i)+s(i)^2-1));


dxx = b(i)/100;

xu = 0:dxx:b(i); xd = b(i)+dxx:dxx:L;

xx = [-flip(xd),-flip(xu),xu,xd];
%xx = [xu,xd];

% Define w-solutions

% % taking into account a
wu = a(i)/(s(i)^2+r(i)-1)+C1(i)*cos(k1*xu)+C2(i)*cos(k2*xu);
wd = a(i)/s(i)^2*(1-exp(-(xd-b(i))))+d(i)*(exp(-s(i)*(xd-b(i)))-exp(-(xd-b(i))));

% without taking into account a

% wu = C1(i)*cos(k1*xu)+C2(i)*cos(k2*xu);
% wd = d(i)*(exp(-s(i)*(xd-b(i)))-exp(-(xd-b(i))));

ww = [flip(wd),flip(wu),wu,wd];

%mass(i) = trapz(xu,wu)+trapz(xd,wd);

% taking into account a

Wu = @(xx) a(i)/(s(i)^2+r(i)-1)+C1(i)*cos(k1*xx)+C2(i)*cos(k2*xx);

Wd = @(xx) a(i)/s(i)^2*(1-exp(-(xx-b(i))))+d(i)*(exp(-s(i)*(xx-b(i)))-exp(-(xx-b(i))));

% without taking into account a 

% Wu = @(xx) C1(i)*cos(k1*xx)+C2(i)*cos(k2*xx);
% 
% Wd = @(xx) d(i)*(exp(-s(i)*(xx-b(i)))-exp(-(xx-b(i))));


mass2(i) = integral(Wu,0,b(i))+integral(Wd,b(i),L);

% Make plot of the solution

% list = [1,100,300,500,700,900,950];
% if ismember(i,list)==1
% figure(i+10)
% plot(xx,ww,'Linewidth',1.5); 
% xlabel('x'); ylabel('w')
% %xlim([xx(1) xx(end)])
% title(['R=',num2str(r(i))])
% end

 end
 


% figure(1)
% set(gcf,'position',[20 50 800 500])


% hfig = figure(1);
% pos = get(hfig,'position');
% set(hfig,'position',pos.*[.5 1 2 1])
figure(1)
x0=10; y0=10; width= 1000; height=400; 
set(gcf,'position',[x0,y0,width,height])

subplot(1,2,1)

plot(R,s_tm,'r','Linewidth',1.5); hold on; 
%plot(R(3:end),s_tm(3:end),'r','Linewidth',1.5); hold on; 
plot(r, s,'b--','Linewidth',1.5); hold on;
xlabel('r'); ylabel('\sigma')
title('\rm Growth Rate')
legend('Finite domain','Infinite domain'); legend boxoff
set(gca,'FontSize',12)
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
annotation('textbox', [0.05, 0.98, 0, 0], 'string', '(a)','fontsize',12)


subplot(1,2,2)

plot(R,b_tm,'r','Linewidth',1.5); hold on;
%plot(R(3:end),b_tm(3:end),'r','Linewidth',1.5); hold on;
plot(r, b,'b--','Linewidth',1.5); hold on;
xlabel('r'); ylabel('b')
title('\rm Half-Ascent Length')
legend('Finite domain','Infinite domain','Location','NorthEast'); legend boxoff
set(gca,'FontSize',12)
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
annotation('textbox', [0.49, 0.99, 0, 0], 'string', '(b)','fontsize',12)

%saveas(gcf,'DRV_Paper/Figures/dispersion','epsc')

% % explore the solution break-down behavior
% figure(2)
% term1 = 1-r.*(2+s.^2);
% term2 = 1-6*r.*s.^2+r.^2.*(s.^4+4*s.^2);
% plot(r,term1); hold on; plot(r,term2); hold on; plot(r,term1-term2)
% xlabel('r'); 
% legend('outer','inner','outer-inner'); legend boxoff
% Plot the residual of the tan-equations

[k1,k2] = Wavenumbers(s,r);

[Eq1,Eq2] = tan_equations(k1,k2,s,b,r);


figure(3)
set(gcf,'position',[20 50 1250 600])
subplot(1,2,1)
loglog(r,abs(Eq1))
xlabel('r');ylabel('residual')
set(gca,'FontSize',14)
title('Tan-Equation1')

subplot(1,2,2)

loglog(r,abs(Eq2))
xlabel('r');ylabel('residual')
set(gca,'FontSize',14)
title('Tan-Equation2')

[k1,k2] = Wavenumbers(s,r);

RHS1 = r.*k1.*k2./(1+s).*(-1./(r.*k2)+s.*k2./(s.^2+r-1));
RHS11 = -k1./(1+s);
RHS12 = r.*s.*k1.*k2.^2./((1+s).*(s.^2+r-1));



RHS2 = r.*k1.*k2./(1+s).*(-1./(r.*k1)+s.*k1./(s.^2+r-1));
RHS21 = -k2./(1+s);
RHS22 = r.*s.*k1.^2.*k2./((1+s).*(s.^2+r-1));


eps= 1e-8;
k0 = 1./sqrt(r).*sqrt(r.^2-3*r+1);
k11 = -0.5*sqrt(r).*(r.^2-3*r+2)./(r.^2-3*r+1);
k21 = 0.5*sqrt(r).*1./(r.^2-3*r+1).^3/2;

k1_approx = k0 + eps*k11;
k2_approx = eps*k21;

figure(100);
plot(r,k1,'b'); hold on;
plot(r,k1_approx,'b--'); hold on;
plot(r,k2,'r'); hold on;
plot(r,k2_approx,'r--');
legend('k1','k_{approx}','k2','k2_{approx}'); legend boxoff