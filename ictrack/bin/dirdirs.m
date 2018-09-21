function dirs = dirdirs(dirname);
% function that returns all the directories under a given directory

dirs = dir(dirname);

isdirs = [dirs.isdir];   % only retain indices of all listings that are directories
isdirs(1) = 0;   % current directory      
isdirs(2) = 0;   % the parent directory 

dirs(~isdirs) = [];
