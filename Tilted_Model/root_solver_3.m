function F = root_solver_3(x,r)

% x(1) = sigma, x(2) = b

k1 = 1/sqrt(2*r)*sqrt(1-r*(x(1)^2+2)+sqrt(1+r^2*(x(1)^4+4*x(1)^2)-6*r*x(1)^2));

k2 = 1/sqrt(2*r)*sqrt(1-r*(x(1)^2+2)-sqrt(1+r^2*(x(1)^4+4*x(1)^2)-6*r*x(1)^2));


% if isreal(k1)==1
% 
% if isreal(k2)== 1 
%  
% [k1,k2] = Wavenumbers(x(1),r);   
%     
% % Full System
% 
% F(1) = tan(k1*x(2))+r*k1*k2/(1+x(1))*(1/(r*k2)-k2*x(1)/(x(1)^2+r-1));
% % 
% F(2) = tan(k2*x(2))-r*k1*k2/(1+x(1))*(-1/(r*k1)+x(1)*k1/(x(1)^2+r-1));
% 
% else 
%   
% [k1,k2] = Wavenumbers_Complex(x(1),r);  
%     
% F(1) = tan(k1*x(2))+r*k1*k2/(1+x(1))*(1/(r*k2)+k2*x(1)/(x(1)^2+r-1));
% % 
% F(2) = tanh(k2*x(2))-r*k1*k2/(1+x(1))*(-1/(r*k1)+x(1)*k1/(x(1)^2+r-1));
% 
% end
% 
% else
%     
% if isreal(k2)==1
%     
% [k1,k2] = Wavenumbers_Complex1(x(1),r);   
%     
% % Full System
% 
% F(1) = tanh(k1*x(2))+r*k1*k2/(1+x(1))*(1/(r*k2)-k2*x(1)/(x(1)^2+r-1));
% % 
% F(2) = tan(k2*x(2))-r*k1*k2/(1+x(1))*(-1/(r*k1)-x(1)*k1/(x(1)^2+r-1));
% 
% else
% [k1,k2] = Wavenumbers_Complex2(x(1),r);   
%     
% % Full System
% 
% F(1) = tanh(k1*x(2))+r*k1*k2/(1+x(1))*(1/(r*k2)+k2*x(1)/(x(1)^2+r-1));
% % 
% F(2) = tanh(k2*x(2))+r*k1*k2/(1+x(1))*(1/(r*k1)+x(1)*k1/(x(1)^2+r-1));
% 
% end
%          
% end
% 
% 

% % % % Solve the two tangent equations
% % % % 
F(1) = tan(k1*x(2))+r*k1*k2/(1+x(1))*(1/(r*k2)-k2*x(1)/(x(1)^2+r-1));
% 
F(2) = tan(k2*x(2))-r*k1*k2/(1+x(1))*(-1/(r*k1)+x(1)*k1/(x(1)^2+r-1));


% if isreal(k2)== 1 
%  
% [k1,k2] = Wavenumbers(x(1),r);   
%     
% % Full System
% 
% F(1) = tan(k1*x(2))+r*k1*k2/(1+x(1))*(1/(r*k2)-k2*x(1)/(x(1)^2+r-1));
% % 
% F(2) = tan(k2*x(2))-r*k1*k2/(1+x(1))*(-1/(r*k1)+x(1)*k1/(x(1)^2+r-1));
% 
% else 
%   
% [k1,k2] = Wavenumbers_Complex(x(1),r);  
%     
% F(1) = tan(k1*x(2))+r*k1*k2/(1+x(1))*(1/(r*k2)+k2*x(1)/(x(1)^2+r-1));
% % 
% F(2) = tanh(k2*x(2))-r*k1*k2/(1+x(1))*(-1/(r*k1)+x(1)*k1/(x(1)^2+r-1));
% 
%    
% end
%     

% Small r approx 1

%F(1) = tan(x(2)/sqrt(r))+1/(sqrt(r)*(x(1)+1));

%F(2) = tan(sqrt(x(1)^2-1)*x(2))-sqrt(x(1)^2-1)*(1-x(1)^2+x(1))/((x(1)+1)*(x(1)^2-1));

% Small r approximation

%F(1) = x(2)-pi/2*sqrt(r)-r*(1+x(1));

%F(2) = x(2)-1/(1+x(1))*(1-x(1)*(x(1)-1))/(x(1)^2-1);


end