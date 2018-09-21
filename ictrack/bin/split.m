function objout = splitobj(objin,X,Y,x1,y1,x2,y2);
% takes in object OBJIN, and coordinates X1,Y1,X2,Y2 that define a line
% then splits OBJ into multiple objects around the dividing line

r = 0:0.01:1;

xp = round(x1 + (x2-x1).*r);   % parametrized version of them line
yp = round(y1 + (y2-y1).*r);   % parametrized version of the line


inds_line = sub2ind([Y X], yp, xp);
inds_obj = sub2ind([Y X], objin.b(:,1), objin.b(:,2));

im = zeros(Y,X);
im(inds_obj) = 1;   % make dividing line;
im = imfill(im,'holes');   % fill holds;
im(inds_line) = 0;

[b,iml] = bwboundaries(im);  % make the boundary image
s = regionprops(iml,'Centroid');

for i = 1:length(s)
    objout(i).m = objin.m;
    objout(i).n = objin.n;
    objout(i).x = s(i).Centroid(1);
    objout(i).y = s(i).Centroid(2);
    objout(i).num = 0;
    objout(i).b = b{i};
    objout(i).trno = 0;    
end


