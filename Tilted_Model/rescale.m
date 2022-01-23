function [LHS] = rescale(Phi)

tol = 10;
ampl = 100;

if max(Phi)>= tol
    
    LHS = Phi/ampl;
else
    LHS = Phi;
    
end

% if rms(Phi)>= tol
%     
%     LHS = Phi/ampl;
%     
% else
%     LHS = Phi;
% end

end