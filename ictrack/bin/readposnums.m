function [M,N] = readposnums(dirname);
% function that reads the number of stage positions in a given micromanager multi-D acquisition dir

pos_expr = 'Pos_(\d\d\d)_(\d\d\d)';  % regular expression for capturing stage positions

% find the total number of stage positions in this sub-directory
posdirs = dir(dirname);
posdirs = [posdirs.name];  % concatenate all the directory names into one large string
posnums = regexp(posdirs, pos_expr, 'tokens');

% first obtain largest x coordinate
xmax = -1;
ymax = -1;

for j = 1:length(posnums)
   xmax = max(xmax, str2num(posnums{j}{1}));
   ymax = max(ymax, str2num(posnums{j}{2}));
end

M = ymax+1;
N = xmax+1;
