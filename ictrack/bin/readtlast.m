function T = readtlast(dname, M, N, channel, z)

C = length(channel);
zslices = length(z.bottom:z.step:z.top);
if (~channel(C).doz);
    lastz = 0;
else
    lastz = zslices-1;
end

%% regular expressions for capturing filenames
img_wild = [dname '/Pos_' num2strn(N-1,3) '_' num2strn(M-1,3) '/img_0000****_' channel(C).name '_' num2strn(lastz,3) '.tif'];
img_expr = ['img_(\d+)_' channel(C).name '_' num2strn(lastz,3) '.tif'];

%% read all the files in the directory
images = dir(img_wild);   % obtain directory listing of all files in the last position
images = [images.name];  % the concatenation of all the files 

%% exit if no files are found
if isempty(images);
    T = 0;
    return
end
ts = regexp(images, img_expr, 'tokens');   % find the last time point
tlast = -1;
for i = 1:length(ts)
    tlast = max(tlast, str2num(ts{i}{1}));
end
T = tlast+1;

