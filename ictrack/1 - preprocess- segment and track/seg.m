function seg(dirname);
% FUNCTION SEG(DIRNAME) takes timelapse images in
% directory DIRNAME, runs segmentation algorithm to identify cells,
% then outputs a file with cell segmentation parameters

%% load the file
load([dirname '/acq.mat']);   % now try get information about the number of z-stacks, the number of stage positions, etc...

%% extract parameters 
p=load('params.mat');
segchannel = p.segchannel;    % 
out = p.segfile;   % the filename for the segmentation file
segfun = p.segfun;  % this is the name of the function used for cell segmentation
clear p;

C = acq.C;
Z = C(segchannel).zslices;
M = acq.M;
N = acq.N;
T = acq.T;
X = acq.X;
Y = acq.Y;
tlist = C(segchannel).tlist;
Xo = 55;  % this is the offset for tiling / tracking
Yo = 55;
objects = [];

%% now process images individually
ind = 0;  % this is the index number for the segmentation frames

for i = 1:T
    clear images;
    if (tlist(i) == 0)
        continue  % skip this frame, no fluroescence data here
    else
        ind = ind+1;
    end
    imfile = [dirname '/imgf_' num2strn(i,4) '.mat'];
    load(imfile);  % loads the images structure
    count = 1;
    
    %% go through all images in the tile
    for m = 1:M
        for n = 1:N                                    
            im = double(images(segchannel).im(:,:,m,n));            
            seg = feval(segfun, im);   % running the cell segmentation algorithm for this given image
            fprintf(['i = ' num2str(i) '\n']);
            [imB, imL] = bwboundaries(seg, 8, 'noholes');     % generates labeled image and boundary list
            s  = regionprops(imL, 'Area','Centroid','PixelIdxList');  % obtain centroids and areas from objects in this one framecentroids and areas from objects in this one frame
            S = length(s);   % total number of objects in this image         
            if (~S)
                objects(ind).obj= []; 
                continue;     % no objects were found, moving to next frame    
            end
            for j = 1:S   
                % make current object
                obj.m = m;
                obj.n = n;
                obj.trno = 0;         
                % segmentation information
                obj.x = s(j).Centroid(1);;
                obj.y = s(j).Centroid(2);
                obj.b = imB{j}; % coordinates for the boundary of the object\                
                % segmentation info to be stored in the data structure
                obj.data.area = s(j).Area;      % object area
                obj.data.perim = length(obj.b);   % object perimeter,     
                obj.data.x = obj.x + (X-Xo).*(obj.n-1);      % real x position
                obj.data.y = obj.y + (Y-Yo).*(obj.m-1);      % real y position with tiling
                obj.gate = 0;    % assign to a null gate                
                % store object into current timepoint and increment the counter
                objects(ind).obj(count) = obj;   % store current object
                count = count+1;
            end
        end
    end
end
gatenames = {'ALL'};
save([dirname '/' out], 'objects', 'gatenames');
