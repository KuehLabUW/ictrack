function [xp, dr, dl] = findxpw(im, drops);

[n,x] = hist(im(:),0:65535);  % extract image histogram
n = n(2:end);  % remove first data point, which contains zeros

% smooth the histogram
smooth = 10;
h = ones(1,smooth)./smooth;
n = filter(h,1,n);

%% find the lo in grayscale intensity
xp = find(n == max(n));
xp = xp(1);  % find the first element


%% now find the hi in grayscale intensity
tail = find(n < max(n)/2);   % find the grayscale value at which intensity drops to half of its maximal value

%% find the location of the right shoulder
tailr = tail;
tailr(tailr < xp) = [];    % remove the tail components that are two the left of xmax
dr = tailr(1)-xp;  % find the first element

%% find the location of the left shoulder (if it exists)

taill = tail;
taill(taill > xp) = [];
dl = xp-taill(end);