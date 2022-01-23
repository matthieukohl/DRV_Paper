% Wavenumbers in the tilted-model
close all;
clear;

N = 1000;
r = 0.01;
sigma = linspace(1,1.5,N);
b = linspace(-0.3,0.3,N);

[S,B] = meshgrid(sigma,b);

k1 = 1/sqrt(2*r)*sqrt(1-r*(S.^2+2)+sqrt(1+r^2*(S.^4+4*S.^2)-6*r*S.^2));

k2 = 1/sqrt(2*r)*sqrt(1-r*(S.^2+2)-sqrt(1+r^2*(S.^4+4*S.^2)-6*r*S.^2));

K1 = k1(11); K2 = k2(11);

k1_approx = 1./sqrt(r)*ones(N,N);
k2_approx = sqrt(S.^2-1);


F2 = tan(k2.*B)-r*k1.*k2./(S+1).*(-1./(r*k1)+k1./S);
F22 = k2.*B-r*k1.*k2./(S+1).*(-1./(r*k1)+k1./S);
%
F1 = tan(k1.*B)+r*k1.*k2./(S+1).*(1./(r*k2)-k2./S);

F11 = tan(k1.*B)+1/sqrt(r)*1./(1+S);

a = 1;

% figure(1)
% contourf(S,B,F1); colorbar
% xlabel('sigma'); ylabel('b')
% title('Equation (*) LHS-RHS')
% % 
% figure(2)
% contourf(S,B,F2); colorbar
% xlabel('sigma'); ylabel('b')
% title('Equation (**) LHS-RHS')
% 
% figure(3)
% plot(B(:,841),F1(:,841),'linewidth',1.6); hold on;
% xlabel('b')
% title(['sigma =',num2str(round(sigma(841),2))])
% legend('Equation (*)')
% legend boxoff
% set(gca,'TickDir','out','Box','off','Layer','top')
% set(gca,'fontsize',14)
% set(gca,'linewidth',1.5)
% 
% figure(4)
% plot(B(:,841),F2(:,841),'linewidth',1.6); hold on;
% xlabel('b')
% title(['sigma =',num2str(round(sigma(841),2))])
% legend('Equation (**)')
% legend boxoff
% set(gca,'TickDir','out','Box','off','Layer','top')
% set(gca,'fontsize',14)
% set(gca,'linewidth',1.5)
% 
% figure(5)
% plot(B(:,810),F1(:,810),'linewidth',1.6); hold on;
% plot(B(:,810),F2(:,810),'linewidth',1.6); hold on;
% xlabel('b')
% title('sigma = 1.40')
% legend('Equation (*)','Equation (**)')
% legend boxoff
% set(gca,'TickDir','out','Box','off','Layer','top')
% set(gca,'fontsize',14)
% set(gca,'linewidth',1.5)
% 
% figure(6)
% plot(B(:,810),F11(:,810),'linewidth',1.6); hold on;
% plot(B(:,810),F22(:,810),'linewidth',1.6); hold on;
% xlabel('b')
% title('sigma = 1.40')
% legend('Equation (*)','Equation (**)')
% legend boxoff
% set(gca,'TickDir','out','Box','off','Layer','top')
% set(gca,'fontsize',14)
% set(gca,'linewidth',1.5)



% 
% figure(3)
% contourf(S,B,F22); colorbar
% xlabel('sigma'); ylabel('b')
% title('tan(k2*b) Approx')
% 
% figure(4)
% contourf(S,B,F11); colorbar
% xlabel('sigma'); ylabel('b')
% title('tan(k1*b) Approx')



% figure(4)
% contourf(S,B,F1-F2); colorbar
% xlabel('sigma'); ylabel('b')
% title('Residual')

% Residual = abs(F1-F2);
% 
% minMatrix = min(Residual(:));
% [row,col] = find(Residual==minMatrix);
% sigma(col)
% b(row)
% Residual(row,col) = NaN;
% minMatrix = min(Residual(:));
% [row,col] = find(Residual==minMatrix);
% sigma(col)
% b(row)

Residual_approx = abs(F11-F22);
minMatrix_approx = min(Residual_approx(:));
[row,col] = find(Residual_approx==minMatrix_approx);
sigma(col)
b(row)
% Residual(row,col) = NaN;
% minMatrix = min(Residual(:));
% [row,col] = find(Residual==minMatrix);
% sigma(col)
% b(row)


% F111 = tan(B/sqrt(r))+1./(sqrt(r)*(1+S));
% F111 = B-pi*sqrt(r)-1./(1+S);
% 
% Res = abs(F111-F22);
% minMatrix = min(Res(:));
% [row,col] = find(Res==minMatrix);
% sigma(col)
% b(row)

