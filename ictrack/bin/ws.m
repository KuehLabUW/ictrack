function objout = ws(objin,X,Y)

objinds = sub2ind([Y X], objin.b(:,1), objin.b(:,2));               
im = zeros(Y,X);   % binary image with the selected objected
im(objinds) = 1;            

im = imfill(im,'holes');   % fill all the holes in the image first
im = -bwdist(~im);
L = watershed(im);
im2 = (L & im);  % this is the watershed segmented image

[b,iml] = bwboundaries(im2);  % make the boundary image
s = regionprops(iml,'Centroid');

for i = 1:length(s)
    objout(i).m = objin.m;
    objout(i).n = objin.n;
    objout(i).x = s(i).Centroid(1);
    objout(i).y = s(i).Centroid(2);    
    objout(i).b = b{i};
    objout(i).trno = 0;  
    objout(i).data = [];
    objout(i).gate = objin.gate;   
end