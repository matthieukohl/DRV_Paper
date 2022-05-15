function [value,isterminal,direction] = myEvent_1(t,y)

tol = 10;

if rms(y(1:length(y)/2))>=tol
    g = 0;
 
else
    g = 1;
end

value = g;  
isterminal = 1;         
direction  = 0;         
end