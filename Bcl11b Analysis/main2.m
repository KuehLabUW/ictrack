% script for running the 
indir = {'/data/phnguyen/Imaging/RawData/082117_p21_EV_cMyc'};
basenames = {'GoodExp'}
outdir = '/data/phnguyen/Imaging/RawData/082117 - processed';

% generate parameter file for image procesing.
paramfile = 'params.mat';
makeparams(paramfile);
pipeline2(indir, basenames, outdir, paramfile);