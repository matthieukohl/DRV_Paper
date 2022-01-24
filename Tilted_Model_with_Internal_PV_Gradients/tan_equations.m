function [RHS1,RHS2] = tan_equations(k1,k2,sigma,b,r)

% Calculate the two tangent equations

RHS1 = tan(k1.*b)+r.*k1.*k2./(1+sigma).*(1./(r.*k2)-sigma.*k2./(sigma.^2+r-1));

RHS2 = tan(k2.*b)-r.*k1.*k2./(1+sigma).*(-1./(r.*k1)+sigma.*k1./(sigma.^2+r-1));

end
