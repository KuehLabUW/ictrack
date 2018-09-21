function [seg] = cellseg122017(im)

log_threshold = -10;  % edge-detection threshold for finding cell boundaries was -10
low_in = 100/600; % lowest intensity input value for cell regonition 
high_in = 300/600; % higest intensity input value to guarantee cell
low_out = 0; % adjusted low intensity output
high_out = 1; % adjusted high intensity output
minsize = 250;   % minimum size of a cell
maxsize = 2000;  % maximum size of a cell
maxecc = 10;    % the maximum eccentricity of a cell was 2
maxratio = 56; % the maximum perimeter squared/area of a cell was 28
% maxwobble = 0.03; % the maximum variation in cell perimeter
unsharp_size = 20;  % 'size' of the unsharp filter was 7
unsharp_alpha = 0.8;  % degree of unsharp filtering between 0.2 and 0.9 was 0.8
% unsharp mask - um - critical for making sharper edges %%
% parameters optimized on 062717 by Blythe Irwin, incorporated negative value for thresholding of LoG image
% changed minimum, maximum cell size, and cell eccentricity
um = fspecial('gaussian',[unsharp_size.*2 unsharp_size.*2],unsharp_size);
h1 = fspecial('gaussian',[7 7], 3); % 2D guassian filter
h2 = fspecial('laplacian', unsharp_alpha);  % what is this shape factor?  worry about this later.
se = strel('disk',2); % morphological structuring element
i1 = imfilter(im,h1);   % gaussian filtered image
i2 = imadjust(i1./600,[low_in; high_in],[low_out; high_out]); % adjust image based on intensity thresholds
i3 = imfilter(600.*i2,h2);   % the laplacian of the gaussian filtered image
i4 = (i3 < -log_threshold);      % select negative values of the laplacian
i5 = imopen(i4,se);   % erodes then dilates image to remove noise
i6 = imclearborder(i5,8);  % clear border objects
i7 = bwmorph(i6,'clean');  % remove isolated pixels
i8 = imfill(i7,'holes');   % fill holes
i9 = bwlabel(i8);   % label objects060117
props = regionprops(i9, im, 'Area','Perimeter','MajorAxisLength','MinorAxisLength','MeanIntensity'); % properties of objects

%figure(n); % display image after each processing step
%subplot(2,3,1); imshow(i1,[]) % original
%subplot(2,3,2); imshow(i2,[]) % gaussian filtered
%subplot(2,3,3); imshow(i3,[]) % laplacian of gaussian
%subplot(2,3,4); imshow(i4,[]) % negative values of laplacian
%subplot(2,3,5); imshow(i5,[]) % eroded and dilated
%subplot(2,3,6); imshow(i6,[]) % clear border objects
%subplot(3,3,7); imshow(i7,[]) % clear border objects
%subplot(3,3,8); imshow(i8,[]) % clear border objects
%subplot(3,3,9); imshow(i9,[]) % clear border objects

% find ways to exclude the non-cells
myecc = [props.MajorAxisLength]./[props.MinorAxisLength]; % eccentricity of each object
area = [props.Area]; % area of each object
ratio = [props.Perimeter].^(2)./[props.Area]; % perimeter squared/area ratio of each object
% [w] = wobble(i9,im); % wobble of each object
select = intersect(find(area > minsize), find(area < maxsize));   % select objects of a certain size
select = intersect(select, find(myecc < maxecc));  % select objects that are not too eccentric
select = intersect(select, find(ratio < maxratio)); % select objects that are below a certain perimeter squared/area ratio 
% select = intersect(select, find(w < maxwobble)); % select objects without too much wobble

seg = ismember(i9, select); % select objects within parameters from labeled objects060117 to output
%[wob] = wobble(seg,im); % wobble of only segmented images
end