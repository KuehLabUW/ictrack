function imout = closest_watershed(im, xnow, ynow)

if (max(im(:)) == 0 )   % empty image
    imout = im;
    return;
end


im = imfill(im,'holes');   % fill all the holes in the image first
im = -bwdist(~im);
L = watershed(im);
im = (L & im);  % this is the watershed segmented image
iml = bwlabel(im);   % this is the binary image


props = regionprops(iml,'Centroid');
N = length(props);   % the number of objects
cs = [props.Centroid];
cs = reshape(cs,2,N)';  % now use this to calculate the object with minimum euclidean distance
sqdist = (cs(:,1) - xnow).^2 + (cs(:,2) - ynow).^2; % calculate the distance matrix from specified coordinate
theone = find(sqdist == min(sqdist));   % identify the index of the one - the closest object 
theone = theone(1);   % only obtain the first index incase there are multiple matches, never know

imout = ismember(iml, theone);  % this is the object
