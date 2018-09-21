function preprocess(indir, basenames, outdir, pos)
% preprocess(BNAMES, INDIR, OUTDIR, I) takes TIF files acquired using
% Metamorph Multi-dimensional acquisition, performs background subtraction
% on the images, and then saves the resultant images as a series of Matlab
% files
% BNAMES gives a cell array of Metamorph base names (in temporal sequence)
% INDIR gives the input directory where TIF files are located 
% OUTDIR gives the output directory
% I represents the index of the position 

% load necessary parameters
p = load('params.mat');
S = p.S;   % stage position

if isstruct(S)  
    mapping = S(pos).mapping;   % this stage position consists of a tile of images
else
    mapping = pos;        % single image, to ensure backwards compatibility
end

M = size(mapping, 1); % the number of row positions
N = size(mapping, 2); % the number of column positions
Z = p.Z; % the number of z sections

% determine structuring element for background subtraction
if isfield(p,'cellsize')
    se = strel('disk',2*p.cellsize);      % flatten with user-defined structuring element
else
    se = strel('disk',50);    % default size of background flattening element
end

segchannel = p.segchannel;  % the fluorescence channel used for segmentation
datachannel = p.datachannel;  % the fluorescence channel used for data skipping 

C = p.C;   % channels structure
% prepare to process all the images
%h = waitbar(0, 'Reading file info...');
B = length(basenames);   % the number of basenames

% find corresponding directories for all the filenames
if (length(indir)==1)&(B>1)
    for i = 2:B
        indir{i} = indir{1};
    end
end

%% loop through all the possible basenames
fprintf(['Retrieving file information for position ' num2str(pos) '...\n']);
for b = 1:B           
    %% obtain listing of all image files in the directory   
    allfiles = dir([indir{b} '/*.TIF']);
    L = length(allfiles);
    imind = 1;
    for i = 1:L     % loop through all files
        [tok, mat] = regexp(allfiles(i).name, [basenames{b} '_w(\d)(.+)_s(\d+)_t(\d+)\.TIF'], 'tokens', 'match');
        if (~isempty(tok))  % this is a true image time point for the metamorph file
            % retrieve all the tokens from the regular expression
            c = str2num(tok{1}{1});   % the channel number
            C(c).name = (tok{1}{2});   % the name of the channel
            s = str2num(tok{1}{3}); % the stage position
            t = str2num(tok{1}{4}); % the time point
            dn = allfiles(i).datenum;   % store the serial date number - information on when the image was acquired            
            name = mat{1};
            bs(b).imfiles(imind).name = name;   % structure containing all filenames
            bs(b).imfiles(imind).c = c;         % channel information
            bs(b).imfiles(imind).s = s;         % stage position
            bs(b).imfiles(imind).t = t;         % time point
            bs(b).imfiles(imind).dn = dn;       % serial date number
            imind = imind + 1;
        end
    end
    %% retrieve and save timestamp information
    ts = [bs(b).imfiles.t];
    bs(b).T = max(ts)-1;   % the number of frames, omit the last one, in case it is incomplete 

    %% obtain information from timestamps
    dns = [bs(b).imfiles.dn];
    for i = 1:bs(b).T    
        tinds = find(ts == i);  % find all indices corresponding to this time point;
        bs(b).tr(i) = min(dns(tinds));    % find the timestamp of the earliest image acquired at this time point;
        % tr stands for time real
    end
    
    if (b == 1)   % choose beginning of the first time point as the elapsed time onset
       tr_min = min([bs(1).tr]);
    end
    bs(b).tr = bs(b).tr - tr_min;
    cs = [bs(b).imfiles.c];
    
    %% find all the time points that exist for a given channel
    for c = 1:length(C)   
        cind = find(cs == c); % find all the indices corresponding to a particular channel    
        texist = sort(unique([bs(b).imfiles(cind).t]));  % sort them and make them unique
        texist(find(texist > bs(b).T)) = [];    %  what is this for?
        tlist = zeros(bs(b).T,1);
        tlist(texist) = 1;
        C(c).tlist = [C(c).tlist ; tlist];
    end
end
tcum = [0 cumsum([bs.T])];   % cumulative sum of all timelapse frames 
%% save acquisition information into structure
acq.M = M;
acq.N = N;
acq.Z = Z;
acq.C = C;
acq.T = tcum(end);
acq.indir = indir;
acq.tr = [bs.tr];   % NEED TO FIX THIS TOMORROW!!!!!
acq.tr_min = tr_min;
acq.segchannel = segchannel;
acq.datachannel = datachannel;

% load single image to obtain image size:
imtemp = readmm(acq, indir{1}, basenames{1},1,1,1,pos);   % load individual image, just to get the size
acq.X = size(imtemp,2);
acq.Y = size(imtemp,1);


%%  save acquisition information into output directory   
mkdir(outdir);
save([outdir '/acq.mat'], 'acq');
    
%% process individual timelapse frames
images = [];
fprintf(['Preprocessing images for pos ' num2str(pos) '...\n']);
for b = 1:B   % loop through different basenames.  If multiple filenames need to be appended           
    for t = 1:bs(b).T   % loop through all the timelapse images
        for c = 1:length(C)   % loop through all channels
            if (C(c).tlist(t+tcum(b)) == 0)  % if frame needs to be skipped
                continue
            else
                for z = 1:C(c).zslices  % loop through all z slices
                    for n = 1:N % loop through all X tiles
                        for m = 1:M % loop through all Y tiles
                             s = mapping(m,n);   % this is the index of the tile location.  important if we want to make a tile
                             im1 = readmm(acq, indir{b}, basenames{b},t,c,z,s);   % load individual image, unprocessed                             
                             if isfield(C,'flatten')   % flatten the image if specified in the settings
                                 if (C(c).flatten)   % flatten the image only if specified
                                     im1 = imtophat(im1,se);
                                 end
                             end
                             
                             if (~isempty(C(c).correct))   % perform illumination correction if subtraction matrix exists...
                                 if isequal(size(im1), size(C(c).correct))   % ... and has correct size
                                     im1 = uint16(double(im1)./C(c).correct);  % background corrected image, converted back to uint16
                                 end
                             end                                                                          
                             images(c).im(:,:,m,n,z) = im1;                             
                        end
                    end
                end
            end
        end
    % save the images structure
    outname = [outdir '/imgf_' num2strn(t+tcum(b),4) '.mat'];
    save(outname, 'images');
    images = [];
end
end
