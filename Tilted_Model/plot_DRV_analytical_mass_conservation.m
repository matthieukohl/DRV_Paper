% Check the mass conservation of the analytical solutions to
% the tan-equations on an infinite domain by
% (1) direct numerical integration of the solutions Int w dx
% (2) integrating by hand the symbolic equations, replacing the aL
% term and then evaluating the terms in the resulting equations
% with progressively larger domains L
% Notes 
% (a) lim L--> inft Int w dx is not the same as inft Int lim L--> w dx
% in this case. Thus the a term needs to be kept in the solutions before
% taking the integrals 
% (b) mass conservation as domain is increased
% (c) Method (2) gives more accurate results (aL term eliminated)

clear; 

L = [32*pi,32*pi*10,32*pi*100,32*pi*1000,32*pi*1e4,32*pi*1e6];

% Define r

r = linspace(1e-3,0.379,1000);

% Define Solver Options
options = optimset('TolX',1e-12,'TolFun',1e-12);
 

% Pre-define values

[s,b,fval,C1,C2,d,a,mass,mass2,mass36,Eq1,Eq2] = deal(zeros(length(r),length(L)));

for j = 1:length(L)
% Define Initial Guess

%x0(1) = 1.5; % first solution
x0(1) = 1.6; % second solution

%x0(2) = 0.2; % first solution
x0(2) = 0.1; % second solution

 for i = 1:length(r)

% linear interpolation as initial guess
  
if i>2
x0(1) = 2*s(i-1,j)-s(i-2,j); 
x0(2) = 2*b(i-1,j)-b(i-2,j); 
end

% Solve the tan-equations for growth rate s and half-ascent area b

fun = @(x)residual_fn(x,r(i));
[x,fval(i)] = fminsearch(fun,x0,options);

s(i,j) = x(1);
b(i,j) = x(2);


% calculate and plot the w-profiles (optional)

% Calculate various coefficients in the solutions

[k1,k2] = Wavenumbers(s(i,j),r(i));

% Calculate the residuals
[Eq1(i,j),Eq2(i,j)] = tan_equations(k1,k2,s(i,j),b(i,j),r(i));

% sign switch for d needed to ensure that w<0 in the descending branch 
% when sigma<1

if s(i,j)>=1
    d(i,j) = 1;
else
    d(i,j) = -1;
end

C1(i,j) = d(i,j)/(cos(k1*b(i,j))*(k2^2-k1^2))*(s(i,j)^2-1)/r(i);
C2(i,j) = d(i,j)/(cos(k2*b(i,j))*(k2^2-k1^2))*(1-s(i,j)^2)/r(i);
a(i,j) = d(i,j)/L(j)*((r(i)-1)*s(i,j)*(s(i,j)-1)/(r(i)+s(i,j)^2-1));


% dxx = b(i,j)/100;
% 
% xu = 0:dxx:b(i,j); xd = b(i,j)+dxx:dxx:L(j);
% 
% xx = [-flip(xd),-flip(xu),xu,xd];
% %xx = [xu,xd];
% 
% %Define w-solutions
% 
% % taking into account a
% wu = a(i,j)/(s(i,j)^2+r(i)-1)+C1(i,j)*cos(k1*xu)+C2(i,j)*cos(k2*xu);
% wd = a(i,j)/s(i,j)^2*(1-exp(-(xd-b(i,j))))+d(i,j)*(exp(-s(i,j)*(xd-b(i,j)))-exp(-(xd-b(i,j))));

% without taking into account a

% wu = C1(i,j)*cos(k1*xu)+C2(i,j)*cos(k2*xu);
% wd = d(i,j)*(exp(-s(i,j)*(xd-b(i,j)))-exp(-(xd-b(i,j))));

% ww = [flip(wd),flip(wu),wu,wd];
% 
% mass(i,j) = trapz(xu,wu)+trapz(xd,wd);

% taking into account a

Wu = @(xx) a(i,j)/(s(i,j)^2+r(i)-1)+C1(i,j)*cos(k1*xx)+C2(i,j)*cos(k2*xx);

Wd = @(xx) a(i,j)/s(i,j)^2*(1-exp(-(xx-b(i,j))))+d(i,j)*(exp(-s(i,j)*(xx-b(i,j)))-exp(-(xx-b(i,j))));

% without taking into account a 

% Wu = @(xx) C1(i,j)*cos(k1*xx)+C2(i,j)*cos(k2*xx);
% 
% Wd = @(xx) d(i,j)*(exp(-s(i,j)*(xx-b(i,j)))-exp(-(xx-b(i,j))));


mass2(i,j) = integral(Wu,0,b(i,j))+integral(Wd,b(i,j),L(j));

% Calculate mass conservation with aL term replaced (eq. 36 paper)
% seems to do better job

mass36(i,j) = a(i,j)*b(i,j)/(s(i,j)^2+r(i)-1)+C1(i,j)/k1*sin(k1*b(i,j))...
    + C2(i,j)/k2*sin(k2*b(i,j))+(r(i)-1)/s(i,j)^2*(a(i,j)*b(i,j)/(s(i,j)^2+r(i)-1)...
    +C1(i,j)/k1*sin(k1*b(i,j))+C2(i,j)/k2*sin(k2*b(i,j)))...
    -a(i,j)/s(i,j)^2+a(i,j)/s(i,j)^2*(exp(-(L(j)-b(i,j)))-1)...
    +d(i,j)*(-1/s(i,j)*(exp(-s(i,j)*(L(j)-b(i,j)))-1)+exp(-(L(j)-b(i,j)))-1);


% Make plot of the solution

% list = [1,100,300,500,750,900,950];
% if ismember(i,list)==1
% figure(i+10)
% plot(xx,ww,'Linewidth',1.5); 
% xlabel('x'); ylabel('w')
% %xlim([xx(1) xx(end)])
% title(['R=',num2str(r(i))])
% end


 end
 
end

% Check mass conservation

figure(1)
set(gcf,'position',[20 50 1250 600])

subplot(1,2,1)
for jj = 1:length(L)
    semilogy(r,abs(mass2(:,jj))); hold on;
end
xlabel('r')
title('Int w')
legend('L','1e1xL','1e2xL','1e3xL','1e4xL','1e6xL','Location','NorthWest')
legend boxoff

subplot(1,2,2)
for jj = 1:length(L)
    semilogy(r,abs(mass36(:,jj))); hold on;
end
xlabel('r')
title('Eq.36')
legend('L','1e1xL','1e2xL','1e3xL','1e4xL','1e6xL','Location','NorthWest')
legend boxoff

 
% Plot the residual of the tan-equations

figure(2)
set(gcf,'position',[20 50 1250 600])
subplot(1,2,1)
loglog(r,abs(Eq1(:,1)))
xlabel('r');ylabel('residual')
set(gca,'FontSize',14)
title('Tan-Equation1')

subplot(1,2,2)

loglog(r,abs(Eq2(:,1)))
xlabel('r');ylabel('residual')
set(gca,'FontSize',14)
title('Tan-Equation2')