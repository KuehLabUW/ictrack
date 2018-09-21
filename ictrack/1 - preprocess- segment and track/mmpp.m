function mmpp(indir, basenames, outdir, offset)
% MMPP(BNAMES, INDIR, OUTDIR, OFFSET) takes TIF files acquired using
% Metamorph Multi-dimensional acquisition, performs background subtraction
% on the images, and then saves the resultant images as a series of Matlab
% files
%
% BNAMES gives a cell array of Metamorph base names (in temporal sequence)
% INDIR gives the input directory where TIF files are located 
% OUTDIR gives the output directory
% OFFSET gives the stage number offset for acquiring other types of data

% first get a listing of all the files in all directories
addpath ../bin

% load necessary parameters
p = load('params.mat');
X = p.X; % the number of row pixels 
Y = p.Y; % the number of column pixels
M = p.M; % the number of row positions
N = p.N; % the number of column positions
Z = p.Z; % the number of z sections
S = p.S;  % the total number of stage positions for this condition.  Do we ever need this?

segchannel = p.segchannel;  % the fluorescence channel used for segmentation
datachannel = p.datachannel;  % the fluorescence channel used for data skipping 

C = p.C;   % channels structure

% prepare to process all the images
%h = waitbar(0, 'Reading file info...');
B = length(basenames);   % the number of basenames

for b = 1:B    
    
    fprintf(['Getting filenames for base number ' num2str(b) '\n']);               
    %% obtain listing of all image files in the directory

    if length(indir)==1
      allfiles = dir([indir{1} '/*.TIF']);
    else
      allfiles = dir([indir{b} '/*.TIF']);
    end

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
            bs(b).imfiles(imind).name = name;             
            bs(b).imfiles(imind).c = c;
            bs(b).imfiles(imind).s = s;
            bs(b).imfiles(imind).t = t;
            bs(b).imfiles(imind).dn = dn;
            imind = imind + 1;
        end
    end
    %% retrieve and save timestamp information
    ts = [bs(b).imfiles.t];
    bs(b).T = max(ts)-1;   % the number of frames, omit the last one, in case it is incomplete 

    %% obtain information on timestamps
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
acq.X = X;
acq.Y = Y;
acq.C = C;
acq.T = tcum(end);
acq.indir = indir;
acq.tr = [bs.tr];   % NEED TO FIX THIS TOMORROW!!!!!
acq.tr_min = tr_min;
acq.segchannel = segchannel;
acq.datachannel = datachannel;

%%  save acquisition information into output directory   
mkdir(outdir);
save([outdir '/acq.mat'], 'acq');
    
%% process timelapse frames
images = [];

for b = 1:B   % loop through different basename
    
    fprintf(['Processing images for base number ' num2str(b) '\n']);
    mapping = offset;   % takes into account possible differing offsets
    
    for t = 1:bs(b).T   % loop through all the timelapse images
        t
        for c = 1:length(C)   % loop through all channels
            if (C(c).tlist(t+tcum(b)) == 0)  % if frame needs to be skipped
                continue
            else
                for z = 1:C(c).zslices  % loop through all z slices
                    for n = 1:N % loop through all X tiles
                        for m = 1:M % loop through all Y tiles
                            
                            if length(indir)==1
                                images(c).im(:,:,m,n,z) = readmm(acq, indir{1}, basenames{b},t,c,z,m,n,mapping);
                            else
                                images(c).im(:,:,m,n,z) = readmm(acq, indir{b}, basenames{b},t,c,z,m,n,mapping);
                            end
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
