function process(indir, basenames, outdir, paramfile)

%% PRE-PROCESS AND TRACK IMAGES %%
%% PROCESSING STARTS
p = load(paramfile);
S = p.S;  % the number of stage positions

for i = 1:S    % loop over all stage positions        
    outdir1 = [outdir '/offset ' num2str(i)];  % output directory for the specific file offset    
    mmpp(indir, basenames, outdir1, i);  % Metamorph preprocessing
    seg(outdir1);    % Cell segmentation
    track(outdir1);    % Cell tracking, based on existing segmentation
end

