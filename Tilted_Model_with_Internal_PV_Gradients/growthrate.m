function [sigma] = growthrate(Phi,t)

% Calculate the growthrate as the average of daily-growthrates over the
% last 6-days of the instability calculation

t_av = 4; t_day = 1;

tt = 1;
ss = 0;

growth = zeros(1,t_av);

for ii = 1:t_av
    
while 1
dt1 = t(end-ss)-t(end-tt);
  
  if dt1>=t_day
      break
  end
  
tt = tt + 1;
end

growth(ii) = 1/dt1*log(rms(Phi(end-ss,:))/rms(Phi(end-tt,:)));

ss = tt;
end

sigma = mean(growth);
end