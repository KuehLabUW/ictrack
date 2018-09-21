function imout = findclosest(im, xnow, ynow)

minotherdist = 3;    % threshold for considering other objects besides the closest one 

if (max(im(:)) == 0 )   % empty image
    imout = im;  % nothing to segment, return an empty image
    return;
end

im = imfill(im,'holes');  % create a filled version of the segmented image
iml = bwlabel(im);   % this is the binary image
props = regionprops(iml,'Centroid');
N = length(props);   % the number of objects
cs = [props.Centroid];
cs = reshape(cs,2,N)';  % now use this to calculate the object with minimum euclidean distance
sqdist = (cs(:,1) - xnow).^2 + (cs(:,2) - ynow).^2; % calculate the distance matrix from specified coordinate

sqdist = sqdist ./ min(sqdist);    % closest object as a normalized distance of unity
theone = find(sqdist == 1);   % identify the index of the closest object - the 'strongest candidate'
theone = theone(1);   % only obtain the first index incase there are multiple matches, never know

% isolate other potential candidates that are of a similar distance away
otherones = setdiff(find(sqdist < minotherdist), theone);    % these are the other candidate matches

imout = ismember(iml, theone) + 2.*ismember(iml,otherones);  % strongest candidate labeled with 1, other potential candidates labeled as 2


