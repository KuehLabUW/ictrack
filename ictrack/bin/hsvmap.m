function map = hsvmap(N, varargin)
% function MAP = HSVMAP(N) returns an MAP, a hsv colormap with N elements

mapall = colormap(hsv);   % HSV colormap
H = length(mapall);       % the number of elements in the colormap
inds = round(1:H/N:H);    % indices for N discrete elements

if (length(varargin) == 0)   % there is an additional argument
    map = mapall(inds,:);
else
    ind = varargin{1};  % particular color index desired
    map = mapall(inds(ind),:);
end