function [value,isterminal,direction] = myEvent(t,y)

tol = 10;

if rms(y)>=tol
    g = 0;
 
else
    g = 1;
end

value = g;  
isterminal = 1;         
direction  = 0;         
end