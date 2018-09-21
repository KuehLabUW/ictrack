function seg = cellseg_edge(im)

log_threshold = 3;  % edge-detection threshold for finding cell boundaries
minsize = 55;   % minimum size of a cell
maxsize = 625;  % maximum size of a cell
maxecc = 5;    % the maximum eccentricity of a cell
unsharp_size = 7;  % 'size' of the unsharp filter
unsharp_alpha = 0.9;  % degree of unsharp filtering between 0.2 and 0.9
% unsharp mask - um - critical for making sharper edges %%
% parameters optimized on 072315 by Sharmayne Siu, incorporated negative value for thresholding of LoG image
% changed minimum, maximum cell size, and cell eccentricity
um = fspecial('gaussian',[unsharp_size.*2 unsharp_size.*2],unsharp_size);
h1 = fspecial('gaussian',[7 7], 3);
h2 = fspecial('laplacian', 0.8);  % what is this shape factor?  worry about this later.
se = strel('disk',1);

i1 = imfilter(im,h1);   % gaussian filtered image
i2 = imfilter(i1,h2);   % the laplacian of the gaussian filtered image
i3 = (i2 < -log_threshold);      % select negative values of the laplacian
i4 = imerode(i3,se);   % erode image to remove noise
i5 = imclearborder(i4);  % clear border objects
i6 = bwmorph(i5,'clean');  % remove isolated pixels
i7 = imfill(i6,'holes');   % fill holes
i8 = bwlabel(i7);   % label objects
props = regionprops(i8, im, 'Area','Perimeter','MajorAxisLength','MinorAxisLength','MeanIntensity');

% find ways to exclude the non-cells
myecc = [props.MajorAxisLength]./[props.MinorAxisLength];
area = [props.Area];
perimeter = [props.Perimeter];
ints = [props.MeanIntensity];
select = intersect(find(area > minsize), find(area < maxsize));   % select objects of a certain size
select = intersect(select, find(myecc < maxecc));  % select images that are not too eccentric

i9 = ismember(i8, select);           % should we dilate the image a little bit?
i10 = imdilate(i9,se);
seg = i10;

end
