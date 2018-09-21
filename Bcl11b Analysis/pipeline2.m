function pipeline2(indir, basenames, outdir, paramfile, varargin)

%% PRE-PROCESS AND TRACK IMAGES %%
%% PROCESSING STARTS
p = load(paramfile);
S = p.S;  % the number of stage positions

% load the optional argument to specify range of stage positions to process
if (~isempty(varargin))
    positions = varargin{1};  % vector listing stage positions to process
else
    positions = 1:S;
end

parfor i = positions    % loop over all stage positions        
    outdir1 = [outdir '/pos ' num2str(i)];  % output directory for the specific file offset    
    %preprocess(indir, basenames, outdir1, i);  % Metamorph preprocessing
    seg(outdir1);    % Cell segmentation
    track(outdir1);  % Cell tracking, based on existing segmentation
end