function [r] = r_factor(w,R)

% Calculate the non-linear r-factor

r = zeros(size(w));

for jj=1:length(w)
    if w(jj)>=0
        r(jj) = R;
    else 
        r(jj) = 1;        
    end
end

end