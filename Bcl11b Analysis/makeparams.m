function makeparams(outfile) 

%% MAKEPARAMS(OUTFILE) generates a parameter file with the filename OUTFILE 
%% for subsequent image conversion, segmentation and tracking for the 
%% immune cell analysis suite.
%% HYK 5/3/2016

%% Image acquisition parameters 
S = 795;    % total number of stage positions
X = 1200;   % the number of row pixels
Y = 1200;   % the number of column pixels
Z = 1;     % the number of z-sections
M = 1;  %sqrt(S); % the number of row tiles
N = 1;   %sqrt(S); % the number of column tiles

segchannel = 2;   % fluorescence channel used for segmentation
datachannel = 2;  % fluorescence channel used to set 'skipping' channel for image viewing

C(1).cR = 255;   % DIC channel parameters
C(1).cG = 255;
C(1).cB = 255;
C(1).doz = 0;
C(1).zslices = 1;
C(1).tlist = [];
C(1).min = 0;
C(1).max = 6000;

C(2).cR = 0;  % Bcl11b-CFP channel parameters
C(2).cG = 255;
C(2).cB = 255;
C(2).doz = 1;
C(2).zslices = 1;
C(2).correct = load('061217_CFP.mat');
C(2).tlist = [];
C(2).min = 30;
C(2).max = 600;

C(3).cR = 0;   % Bcl11b-YFP channel parameters
C(3).cG = 255;
C(3).cB = 0;
C(3).doz = 1;  
C(3).zslices = 1;
C(3).correct = load('061217_YFP.mat');   
C(3).tlist = [];
C(3).min = 30;
C(3).max = 400;

C(4).cR = 255;   % Bcl11b-mCherry channel parameters
C(4).cG = 0;
C(4).cB = 0;
C(4).doz = 1;  
C(4).zslices = 1;
C(4).correct = load('061217_mCherry.mat');   
C(4).tlist = [];
C(4).min = 30;
C(4).max = 400;


%% MATLAB function for cell segmentation
segfun = 'cellseg122017';  % user-written function for cell segmentation

%% Object tracking parameters
Cxy = 20;    % maximum travel distance
Carea = 0.6;   % maximal fold change in area
Cf = 0.5;     % maximal fold change in red or green fluorescence level
maxskip = 5; % the total possible number of omitted frames

%% outfile filenames:
segfile = 'segment.mat';    % name for MATLAB data file containing all the segmentation data
segtrackfile = 'segtrack.mat';   % name for MATLAB file contaning segmentation and tracking data

save(outfile);