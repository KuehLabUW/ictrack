function [xnow, ynow] = gcaxy(h);
currentxy = get(h,'CurrentPoint');
xnow = round(currentxy(1,1));
ynow = round(currentxy(1,2));

end

