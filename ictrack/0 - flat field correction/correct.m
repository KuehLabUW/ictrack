function correct(indir, base, channel, out)

%%%% flat field correction for metamorph generated 12 x 12 grid of bead
%%%% positions %%%%% 

%%% updated 08/24/10 %%%%%%% 

%% clear existing variables
clear imname;
clear imfiles;
addpath ../bin

%% multi-dimensional acquisition information for the files
%indir= '/home/kueh/012913 pu1gfp with 021513';
%base = 'beads1';   % the base name for the files
%channel = 'w3GFP';   % the name of the channel
%out = '071811GFP.mat'

N = 100;

%% image segmentation parameters
se = strel('disk',40);
X = 512;
Y = 512;
ints = [];
yall = [];
xall = [];

imall = zeros(Y,X);
imcount = zeros(Y,X);

%% load and segment all the images
fprintf('loading and segmenting images...');
for i = 1:N     
        %% load and segment image
        filename = [indir '/' base '_' channel '_s' num2str(i) '.TIF'];
        if isempty(dir(filename))   % file not found
            continue
        end
        
        im = imread(filename);               
        im(im > 63000) = 0;   % remove saturated pixels 
     
        im = imtophat(im,se);   % background subtraction 
        level = graythresh(im); % compute threshold using Otsu's method
        imb = im2bw(im,level);    % threshold image
        
        imall = imall + (double(imb).*double(im));
        imcount = imcount + imb;                    
end

imcount(find(~imcount)) = NaN;   % x and y positions without data points will not be counted
imavg = imall./imcount;    % here is the averaged x and y intensities at all positions

iminds = find(~isnan(imcount));

[yall, xall] = ind2sub([Y,X], iminds);
intsall = imavg(iminds);

% Convert all inputs to column vectors.
xall = xall(:);
yall = yall(:);
intsall = intsall(:);

ft = fittype( 'poly22' );
opts = fitoptions( ft );
opts.Lower = [-Inf -Inf -Inf -Inf -Inf -Inf];
opts.Robust = 'Bisquare';
opts.Upper = [Inf Inf Inf Inf Inf Inf];
opts.Weights = zeros(1,0);
[fitresult, gof] = fit( [xall, yall], intsall, ft, opts );

[x,y] = meshgrid(1:X, 1:Y);
z = fitresult(x,y);   % evalute the function
z2 = z./max(z(:))  % normalize so that the maximum of the fit is unity
%z2 = z2(2:2:end, 2:2:end);    % correct for difference in binning, not
%required any more
surf(z2); shading interp;
save(out, 'z2');

