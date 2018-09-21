function imout = readmm(acq, indir, base,t,c,z,s);

%indir = acq.indir;   % name of the input directory of the images
cname = ['w' num2str(c) acq.C(c).name];   % name of the channel in Metamorph format 
correct = acq.C(c).correct;

% reads an individual metamorph file from a multi-dimensional
% acquisition

filename = [base '_' cname '_s' num2str(s) '_t' num2str(t) '.TIF'];  % this is the file!
imout = imread([indir '/' filename]);   % load raw image

% do the tophat only for the fluorescence channels
% if (c > 1)
%     se = strel('disk',25);      % flatten with 
%     imout = imtophat(im, se);   % tophat flattening of background
% else
%     imout = im;
% end
% 
% if (~isempty(correct))   % background flattening
%     correct = correct.z2;        
%     imout = uint16(double(imout)./ correct);  % background corrected version of image
% end
% end

