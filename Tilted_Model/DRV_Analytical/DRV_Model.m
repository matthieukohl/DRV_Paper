% Wavenumbers in the tilted-model
close all;
clear;

N = 1000;
r = linspace(0.001,1,N);
R = 0.01;
sigma = 1.41;
b = linspace(-0.3,0.3,N);

k1 = 1./sqrt(2*r).*sqrt(1-r*(sigma^2+2)+sqrt(1+r.^2*(sigma^4+4*sigma^2)-6*r*sigma^2));

k2 = 1./sqrt(2*r).*sqrt(1-r*(sigma^2+2)-sqrt(1+r.^2*(sigma^4+4*sigma^2)-6*r*sigma^2));

K1 = k1(11); K2 = k2(11);

k1_approx = 1./sqrt(r);
k2_approx = sqrt(sigma^2-1)*ones(1,N);%+(4*sigma^4-sigma^2/2)*r/(4*(sigma^2-1));


%growth = (-pi*r-2*sqrt(r)+sqrt(r.*(pi^2*r+12*pi*sqrt(r)+4)))./(2*pi*r);

%plot(r,growth)



F1 = tan(K2*b)-R*K1*K2/(sigma+1)*(-1/(R*K1)+K1/sigma);
F11 = K2*b-R*K1*K2/(sigma+1)*(-1/(R*K1)+K1/sigma);

F2 = tan(K1*b)+R*K1*K2/(sigma+1)*(1/(R*K2)-K2/sigma);

figure(20)
plot(b,F1); hold on;
plot(b,F2); 
xlabel('b')
legend('tan(k2)','tan(k1)')


figure(1)
semilogx(r,k1); hold on;
semilogx(r,k1_approx);
xlabel('r')
ylabel('k1')
legend('k1','k1_{approx}')

figure(2)
semilogx(r,k2); hold on;
semilogx(r,k2_approx)
legend('k2','k2_{approx}')


% x = linspace(b,10,100);
% w_descent = exp(-sigma*(x-b))-exp(-(x-b));
% 
% figure(3)
% plot(x,w_descent)