function imout = readmmorph(acq,base,t,c,z,m,n,mapping);

indir = acq.indir;   % name of the input directory of the images
cname = ['w' num2str(c) acq.C(c).name];   % name of the channel in Metamorph format 
correct = acq.C(c).correct;

% reads an individual metamorph file from a multi-dimensional
% acquisition

% tweak the mapping 082410 to handle side by side conditions
%mapping = [7 8 9 16 17 18 ; 4 5 6 13 14 15;1 2 3 10 11 12];   % for a 3x3
%grid
%mapping = [13 14 15 16 29 30 31 32; 9 10 11 12 25 26 27 28 ;  5 6 7 8 21 22 23 24; 1 2 3 4 17 18 19 20];

s = mapping(m,n);
filename = [base '_' cname '_s' num2str(s) '_t' num2str(t) '.TIF'];  % this is the file!
im = imread([indir '/' filename]);   % load this one particular z plane
% disabled multi-page tiff reading 10/29/2010
% now process the file by flattening it!
se = strel('disk',20); 
imout = imtophat(im, se);   % tophat flattening of background

if (~isempty(correct))
    correct = correct.z2;
    imout = uint16(double(imout)./ correct);  % background corrected version of image
end
end

