function y = residual_fn_without_complex_penalty(x,r)

s = x(1);
b = x(2);

[k1,k2] = Wavenumbers(s,r);

% Define the tan-equations 
y1 = tan(k1*b)+r*k1*k2/(s+1)*(1/(r*k2)-s*k2/(s^2+r-1));
y2 = tan(k2*b)-r*k1*k2/(s+1)*(-1/(r*k1)+s*k1/(s^2+r-1));
 
y = real(abs(y1).^2+abs(y2)^2); 