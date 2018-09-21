function data = cellprops(images, obj)   
% FUNCTION DATA = CELLPROPS(IMAGES, OBJ) takes an input an ictrack image
% matrix IMAGE, an segmented objected OBJ, and returns a DATA structure 
% containing relevant experiment-specific information about the properties
% of the cell object.
%
% This CELLPROPS function is relevant for measuring nuclear/cytoplasmic
% ratios of Erk-BFP (channel 2), H2B-IFP (channel 3) and NFAT-Ruby
% (channel 4).  IFP and Ruby are used for segmenting nuclei and cytoplasm
% respectively.  Canny edge detection with minimal contrast saturation used
% to accurate discern nuclear edges
% Modified 061718, HYK

bfp = images(2).im;   % bfp image
gfp = images(3).im;   % gfp image
ifp = images(5).im;   % H2B-IFP
rfp = images(4).im;   % nfat-ruby

[X,Y] = size(ifp);    % maximum size of the image

% parameters and filters
bbscaling = 2.5;    % enlargement factor for bounding box
se = strel('disk',3);   % structuring element for filling in holes
f1 = fspecial('gaussian',[3 3], 1); 

% first let's display the nucleus, along with the segmentation
% figure(5); subplot(1,3,1); imshow(ifp,[]); hold on; plot(obj.b(:,2), obj.b(:,1),'r-'); title('h2b-ifp')% old segmentation
% subplot(1,3,2); imshow(bfp,[]); hold on; plot(obj.b(:,2), obj.b(:,1),'r-'); title('erkktr-bfp');
% subplot(1,3,3); imshow(rfp,[]); hold on; plot(obj.b(:,2), obj.b(:,1),'r-'); title('nfat-ruby');

% get and plot bounding box for the old segmentation
x1 = min(obj.b(:,2));
x2 = max(obj.b(:,2));
y1 = min(obj.b(:,1));
y2 = max(obj.b(:,1));

for i = 1:3
    subplot(1,3,i);    
    plot([x1 x1 x2 x2 x1],[y1 y2 y2 y1 y1],'r:');
end

xm = (x1+x2)/2;   % mean x position
ym = (y1+y2)/2;   % mean y position

xmin = max(round(xm-(xm-x1)*bbscaling),1);
ymin = max(round(ym-(ym-y1)*bbscaling),1);
xmax = min(round(xm+(x2-xm)*bbscaling),X);
ymax = min(round(ym+(y2-ym)*bbscaling),Y);

% for i = 1:3
%     subplot(1,3,i);    
%     plot([xmin xmin xmax xmax xmin],[ymin ymax ymax ymin ymin],'g:');
% end

% now let us take region of the original image
ifp2 = ifp(ymin:ymax, xmin:xmax);
gfp2 = gfp(ymin:ymax, xmin:xmax);
rfp2 = rfp(ymin:ymax, xmin:xmax);
bfp2 = bfp(ymin:ymax, xmin:xmax);

% figure(6); subplot(1,3,1); imshow(ifp2,[]);
% subplot(1,3,2); imshow(bfp2,[]);
% subplot(1,3,3); imshow(rfp2,[]);

%% Segment H2B-IFP using canny's method
ifp3 = mat2gray(imfilter(ifp2,f1));  % do a gaussian filter, stretch to between 0 and 1     
imin = quantile(ifp3(:),0.75);
imax = max(ifp3(:));
ifp4 = imadjust(ifp3,[imin imax],[0 1]);
ifp5 = edge(ifp4,'canny',0.2);
ifp6 = imclose(ifp5, se); 
ifp7 = imfill(ifp6,'holes');
ifp8 = imopen(ifp7, se);

% plot segmentation
% figure(7); subplot(2,3,1); imshow(ifp3,[]);  title('h2b-ifp');
% subplot(2,3,2); imshow(ifp4,[]);
% subplot(2,3,3); imshow(ifp5,[]);
% subplot(2,3,4); imshow(ifp6,[]);
% subplot(2,3,5); imshow(ifp7,[]);
% subplot(2,3,6); imshow(ifp8,[]);

%% Segment NFAT-Ruby using canny's method
rfp3 = mat2gray(imfilter(rfp2,f1));  % do a gaussian filter, stretch to between 0 and 1     
imin = quantile(rfp3(:),0.75);
imax = max(rfp3(:));
rfp4 = imadjust(rfp3,[imin imax],[0 1]);
rfp5 = edge(rfp4,'canny',0.2);
rfp6 = imclose(rfp5, se); 
rfp7 = imfill(rfp6,'holes');
rfp8 = imopen(rfp7, se);

% plot segmentation
% figure(8); subplot(2,3,1); imshow(rfp3,[]);   title('nfat-ruby')
% subplot(2,3,2); imshow(rfp4,[]);
% subplot(2,3,3); imshow(rfp5,[]);
% subplot(2,3,4); imshow(rfp6,[]);
% subplot(2,3,5); imshow(rfp7,[]);
% subplot(2,3,6); imshow(rfp8,[]);

%% Define and plot the nuclear and cytoplasmic masks
nbw = ifp8;
cbw = rfp8-ifp8;
% figure(10); subplot(1,2,1); imshow(nbw,[]); title('nuclear');
% subplot(1,2,2); imshow(cbw,[]); title('cytoplasmic');

% now, we need to save all the data - use the median to avoid outliers due
% outliers due to aberrant values at the boundary.
data.bfp_n = double(median(bfp2(find(nbw))));
data.bfp_c = double(median(bfp2(find(cbw))));
data.gfp_n = double(median(gfp2(find(nbw))));
data.gfp_c = double(median(gfp2(find(cbw))));
data.ifp_n = double(median(ifp2(find(nbw))));
data.ifp_c = double(median(ifp2(find(cbw))));
data.rfp_n = double(median(rfp2(find(nbw))));
data.rfp_c = double(median(rfp2(find(cbw))));
data.bfp_r = data.bfp_n / (data.bfp_n + data.bfp_c);
data.ifp_r = data.ifp_n / (data.ifp_n + data.ifp_c);
data.rfp_r = data.rfp_n / (data.rfp_n + data.rfp_c);
data.gfp_r = data.gfp_n / (data.gfp_n + data.gfp_c);
data.nbw = nbw;
data.cbw = cbw;
data.ifp2 = ifp2;
data.rfp2 = rfp2;
data.bfp2 = bfp2;









