% Curve fitting w to w= mean(r.*w)/s^2+c1*cos(k1*(x-x_max))+c2*cos(k2*(x-x_max)) gives:

close all; clear;

% Load Time-marching solution

load('DRV.mat','x','dx','N','s','r',...
    'w','phi_end','tau_end','partial_t_term')

w = w/norm(w); R = 0.01;

k1 = 1/sqrt(2*R)*sqrt(1-R*(2+s^2)+sqrt(1+R^2*(s^4+4*s^2)-6*R*s^2));

k2 = 1/sqrt(2*R)*sqrt(1-R*(2+s^2)-sqrt(1+R^2*(s^4+4*s^2)-6*R*s^2));


index = find(w>0); 
%Ascent region
xx = x(index(1):index(end)); yy = w(index(1):index(end));

%Descent Region
xx = x(index(end)+1:end-1)'; yy = w(index(end)+1:end);



center = xx(1)+(xx(end)-xx(1))/2;

c1 = 0.3061; c2 = 0.07002; k1 = 9.7727; k2 = 0.8771; x_max =  7.1251;

constant = -0.0019; b = 0.1885; bb = 7.3136; d2 = 0.08629; 

a = -0.0034; r = 0.01; s = 1.3133; L = 8*pi;

% Boundary Equations: Finite Domain

F(1) = a/(r+s^2-1)+c1*cos(k1*b)+c2*cos(k2*b);

F(2) = -r*c1*k1*sin(k1*b)-r*c2*k2*sin(k2*b)-a/s^2-d2*(1-s);

F(3) = a*b/(r+s^2-1)+c1/k1*sin(k1*b)+c2/k2*sin(k2*b)+...
    a/s^2*(L-b)+a/s^2*(exp(-(L-b))-1)+...
    d2*(-1/s*(exp(-(L-b))-1)+(exp(-(L-b))-1));

F(4) = -r*c1*k1^2*cos(k1*b)-r*c2*k2^2*cos(k2*b)+a/s^2-d2*(s^2-1);

F(5) = a*b*r/(r+s^2-1)+r*c1/k1*sin(k1*b)+r*c2/k2*sin(k2*b)+...
    a/s^2*(L-b)+a/s^2*(exp(-(L-b))-1)+...
    d2*(-1/s*(exp(-(L-b))-1)+(exp(-(L-b))-1))-a*L;

% Boundary Equations: Infinite Domain

G(1) = c1*cos(k1*b)+c2*cos(k2*b);

G(2) = -r*c1*k1*sin(k1*b)-r*c2*k2*sin(k2*b)-d2*(1-s);

G(3) = c1/k1*sin(k1*b)+c2/k2*sin(k2*b)-d2*s*(s-1)/(s^2+r-1);

G(4) = -r*c1*k1^2*cos(k1*b)-r*c2*k2^2*cos(k2*b)-d2*(s^2-1);

% Comparison to DRV-Theory: Finite Domain

L = N*dx; r = 0.01;
x0 = [1.3133,0.1885,0.3061,0.07002,-0.0034]; 
sol = fsolve(@(x) root_solver(r,L,x),x0);
sigma = sol(1); b = sol(2); c1 = sol(3); c2 = sol(4); a = sol(5);

k1 = 1/sqrt(2*r)*sqrt(1-r*(2+sigma^2)+sqrt(1+r^2*(sigma^4+4*sigma^2)-6*r*sigma^2));

k2 = 1/sqrt(2*r)*sqrt(1-r*(2+sigma^2)-sqrt(1+r^2*(sigma^4+4*sigma^2)-6*r*sigma^2));

index = find(w>0); 

xu = linspace(-sol(2),sol(2),length(index)+2); 
xd = linspace(sol(2)+dx,sol(2)+dx+250*dx,252); 


wu = a/(r+sigma^2-1)+c1*cos(k1*xu)+c2*cos(k2*xu);
wd = a/sigma^2*(1-exp(-(xd-b)))+(exp(-sigma*(xd-b))-exp(-(xd-b)));

xx = [-flip(xd),xu,xd]; ww = [flip(wd),wu,wd];

%% Finite Domain
% figure(1)
% plot(xx,w(index(1)-1-252:index(end)+1+252),'Color','red','linewidth',1.6); hold on;
% plot(xx,ww/norm(ww),'--','Color','blue','linewidth',1.6); 
% xlabel('x');
% ylabel('w');
% legend('Numerical','Theory')
% set(gca,'fontsize', 12);
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% legend boxoff
% xlim([-max(xd),max(xd)])
% %saveas(gcf,'DRV','epsc')
% 
% Comparison to DRV-Theory: Infinite Domain

r = 0.01;

x0 = [1.3133,0.1885,0.3061,0.07002];
sol = fsolve(@(x) root_solver_2(x,r),x0);

%x00 = [1.3133,0.1885];
%sol = fsolve(@(x) root_solver_3(x,r),x00);

sigma = sol(1); b = sol(2); c1 = sol(3); c2 = sol(4);
index = find(w>0); 

k1 = 1/sqrt(2*r)*sqrt(1-r*(2+sigma^2)+sqrt(1+r^2*(sigma^4+4*sigma^2)-6*r*sigma^2));

k2 = 1/sqrt(2*r)*sqrt(1-r*(2+sigma^2)-sqrt(1+r^2*(sigma^4+4*sigma^2)-6*r*sigma^2));


xu = linspace(-sol(2),sol(2),length(index)+2); 
xd = linspace(sol(2)+dx,sol(2)+dx+250*dx,252); 


wu = c1*cos(k1*xu)+c2*cos(k2*xu);
wd = (exp(-sigma*(xd-b))-exp(-(xd-b)));

xx = [-flip(xd),xu,xd]; ww = [flip(wd),wu,wd];
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


% Check Asymptotic Scaling

M = 100;
r = linspace(0.001,0.1,M);
sigma = zeros(M,1); b = zeros(M,1);

% Theory: Parameters
c0 = 1/2*(1+sqrt(5)); c1 = pi/2*(1+c0)*(c0^2-1)/(1-2*c0);
d2 = c1+pi*(1+2*c0^2)/4;

lambda = c0^2*(2*c0^2-1)/(2*sqrt(c0^2-1));

RHS = c0^2*(1-2*c0^2)-c0*(1+2*c0^2)-c0/(c0^2-1)+c0^3*(2*c0^2-1)/(2*(c0^2-1))+...
    lambda-lambda*c0/(c0^2-1);


c2 = (pi/2*(c0^2*c1+2*c0*c1)+(1+c0)^2*(c0^2-1)+c1^2)/(1-2*c0); %(most likely incorrect
% but for some reason including this term is useful to have better initial
% guess ...)

%c2 = 1;

%s_theory = 1.6180-2.9756*sqrt(r);
%b_theory = pi/2*sqrt(r)+(1+1.6180)*r;

s_theory = c0+c1*sqrt(r);
b_theory = pi/2*sqrt(r)+(1+c0)*r+d2*r.^(3/2);

% s_theory = c0+c1*sqrt(r)+c2*r;
% b_theory = pi/2*sqrt(r)+(1+c0)*r;

for ii = 1:M
x0 = [c0+c1*sqrt(r(ii))+c2*r(ii);pi/2*sqrt(r(ii))+(1+c0)*r(ii)+d2*r(ii).^(3/2)];
%x0 = [rand,rand];
sol = fsolve(@(x) root_solver_3(x,r(ii)),x0);
sigma(ii) = sol(1); b(ii) = sol(2);
end

% Time Marching Approach

r_time = [0.1000,0.0800,0.0600,0.0400,0.0100,0.0080,0.0060,0.0040,0.0010];
b_time = [0.9048,0.7791,0.6283,0.4524,0.2011,0.1759,0.1508,0.1257,0.0503];
s_time = [1.0419,1.0811,1.1290,1.1914,1.3620,1.3836,1.4092,1.4414,1.5247];

close all;
figure(1)
semilogx(r_time,s_time); hold on;
semilogx(r,sigma); hold on;
semilogx(r,s_theory); hold on;
xlabel('r'); ylabel('sigma'); legend('Time Marching','Numerical','Theory');
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);
set(gca,'linewidth',1.5)
%saveas(gcf,'sigma','epsc')
%close all;


figure(2)
semilogx(r_time,b_time); hold on;
semilogx(r,b); hold on;
%semilogx(r,b_theory); hold on;
xlabel('r'); ylabel('b'); legend('Time Marching','Numerical','Theory');
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);
set(gca,'linewidth',1.5)
%saveas(gcf,'b','epsc')

% Test various approximations
N = 100;
r = linspace(0.001,0.1,N); s = 1.3;


k1 = 1./sqrt(2*r).*sqrt(1-r*(2+s^2)+sqrt(1+r.^2*(s^4+4*s^2)-6*r*s^2));
 
k2 = 1./sqrt(2*r).*sqrt(1-r*(2+s^2)-sqrt(1+r.^2*(s^4+4*s^2)-6*r*s^2));

k1_approx = 1./(sqrt(r))-(1+2*s^2)*sqrt(r)/2;
k2_approx = sqrt(s^2-1)+s^2*(2*s^2-1)*r/(2*sqrt(s^2-1))+...
    s^4*(20*s^4-32*s^2+11)*r.^2/(8*(s^2-1)^(3/2));

close all;
figure(1)
plot(r,k1); hold on; plot(r,k1_approx)


figure(2)
plot(r,k2); hold on; plot(r,k2_approx)



% xcentre = x(index(1))+(x(index(end))-x(index(1)))/2;
% 
% load('DRV_40pi.mat','x','dx','N','s','r',...
%     'w','phi_end','tau_end','partial_t_term')

%load('DRV_80pi.mat','x','dx','N','s','r',...
%'w','phi_end','tau_end','partial_t_term')

% load('DRV_80pi.mat','x','dx','N','s','r',...
% 'w','phi_end','tau_end','partial_t_term')

% load('DRV_16pi.mat','x','dx','N','s','r',...
%     'w','phi_end','tau_end','partial_t_term')
% 
% index = find(w>0); bb = (x(index(end)+1)-x(index(1)-1))/2;

