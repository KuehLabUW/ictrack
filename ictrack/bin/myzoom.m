function myzoom(factor,maxzoom);

% POS - current cursor position
% FACTOR - the degree of zoom
% IN - the current axis scaling
% LIM - the scaling of the maximally zoomed out image

in = [get(gca,'XLim') get(gca,'YLim')];
lim = [0 maxzoom 0 maxzoom];
[x,y] = gcaxy(gca);
pos = [x y];

dX_max = lim(2) - lim(1);
dY_max = lim(4) - lim(3);

dX_in = in(2) - in(1);
dY_in = in(4) - in(3);

dX_out = min(dX_max, dX_in / factor);
dY_out = min(dY_max, dY_in / factor);

out(1) = pos(1) - 0.5*dX_out;
out(2) = pos(1) + 0.5*dX_out;
out(3) = pos(2) - 0.5*dY_out;
out(4) = pos(2) + 0.5*dY_out;

if out(1) < 0
    out(1) = 0;
    out(2) = dX_out;
elseif out(2) > dX_max
    out(2) = dX_max;
    out(1) = dX_max - dX_out;
end

if out(3) < 0
    out(3) = 0;
    out(4) = dY_out;
elseif out(4) > dY_max
    out(4) = dY_max;
    out(3) = dY_max - dY_out;
end
axis(out);





