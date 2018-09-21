function [xlow, xhigh] = imadjustf(im,bg_drop)

[n,x] = hist(im(:),0:max(im(:)));  % extract image histogram
n = n(2:end); x = x(2:end);  % remove first data point, which contains zeros

%% based on image histogram, adjust contrast in the image
xlow = x(find(n == max(n)));
xlow = xlow(1);  % find the first element    % will use this value later!!!!!!!
outliers = find(n < bg_drop*n(1));
xo = x(outliers);    % the intensity values
no = n(outliers);
xhigh = sum(xo.*no)/sum(no);

return