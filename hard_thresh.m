% Hard Thresholding
function [v]= hard_thresh(x,beta_thresh)
for i = 1:1:length(x)
 if abs(x(i) <= beta_thresh)
 x(i) = 0.0;
 end
 v(i) = x(i);
end
end
