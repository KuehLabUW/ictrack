function timepoints = track(dirname)
% TIMEPOINTS = MUNKRES(IN, OUT)  returns a
% timepoints objects structure with tracking numbers corresponding to
% individual lineages.
% CXY   - maximum travel distance
% CAREA - maximal fold change in area
% CGR   - maximal fold change in red or green fluorescence level
% MAXSKIP - maximal number of ommitted frames
% created, HYK 6/24/11


% Recursion limit for the Hungarian Algorithm
set(0,'RecursionLimit',50000);

% define input and output parameter filenames
p = load('params.mat');
in = [dirname '/' p.segfile];    % the filename for the segmentation data only
out = [dirname '/' p.segtrackfile];  % the filename for the tracking data

% load tracking parameters 
Cxy = p.Cxy;       % maximum change in distance
Carea = p.Carea;    % maximal fold change in area
Cf = p.Cf;   %0.5;     % maximal fold change in fluorescence level
maxskip = p.maxskip;   %1; % the total possible number of ommitted frames

% load the segmentation files themselves
load(in);      % this is the MATLAB data file containing the segmented objects
T = length(objects);   % obtain total number of timepoints from data

%% loop through all pairs of tracks
Ls = [];   % lengths of all the tracks;
for i = 1:(T-1)
    i
    obj1 = objects(i).obj;       % object list for current frame
    obj2 = objects(i+1).obj;     % object list for subsequent frame
    
    if (i == 1)   % initialize objnow structure for tracking, if objects are present
        if length(obj1) > 0
            for m = 1:length(obj1)
                objnow(m).t = i;        % this is the timepoint number, useful later for omitted tracks
                objnow(m).ind = m;      % this is the index number
                objnow(m).tr = m;       % track number
                objnow(m).skipped = 0;
                Ls = [Ls 1];
            end
        else
            objnow = [];
        end
    else    % create new objnow structure for this frame, carrying over all untracked objects
        if exist('objcarry');
            objnow = objcarry;      % these are the untracked objects from prior frames
        else
            objnow = [];
        end
        
        for m = 1:length(obj1);
            if (~obj1(m).trno)   % object has not yet been assigned to a track, create a new track
                Ls = [Ls 1];     % this is the vector of track lengths
                obj1(m).trno = length(Ls);     % create a new track
            end
            objnow(end+1).t = i;   % append new object to the end
            objnow(end).ind = m;
            objnow(end).tr = obj1(m).trno;
            objnow(end).skipped = 0;
        end
    end
  
    %% calculate cost matrix
    % compute cost matrix
    M = length(objnow);  %munkres.m
    N = length(obj2);
    
    % these are all the variables that go into the cost function
    xs1 = zeros(M,N);
    ys1 = zeros(M,N);
    areas1 = zeros(M,N);
    gates1 = zeros(M,N);
    
    xs2 = zeros(M,N);
    ys2 = zeros(M,N);
    areas2 = zeros(M,N);
    gates2 = zeros(M,N);
    
    for m = 1:M      % compute inputs into cost matrix for initial frame
        t1 = objnow(m).t;   % this is the timepoint and the index information for the current timepoint
        ind1 = objnow(m).ind;
        xs1(m,:) = objects(t1).obj(ind1).data.x;
        ys1(m,:) = objects(t1).obj(ind1).data.y;
        areas1(m,:) = objects(t1).obj(ind1).data.area;
        gate = objects(t1).obj(ind1).gate;
        gates1(m,:) = gate;
    end
    
    t2 = i+1;
    for n = 1:N        % cost inputs into cost matrix for subsequent frame
        ind2 = n;
        xs2(:,n) = objects(t2).obj(ind2).data.x;
        ys2(:,n) = objects(t2).obj(ind2).data.y;
        areas2(:,n) = objects(t2).obj(ind2).data.area;
        gate = objects(t2).obj(ind2).gate;
        gates2(:,n) = gate;
        
    end
    
    Dxy = sqrt((xs2-xs1).^2 + (ys2-ys1).^2)./Cxy;  % calculate cost function
    Dxy(find(Dxy > 1)) = Inf;
    Dareas = abs(log2(areas2./areas1));   % modified 072715 - changed expression to measure fold change in area
    Dareas(find(Dareas > Carea)) = Inf;    
    Dgates = abs(gates2-gates1);
    Dgates(find(Dgates > 0)) = Inf;   % forbidden assignment for anything that does not fall into the same gate    
    cost = Dxy + Dareas + Dgates;
    
    %% perform assignment using cost matrix
    assign = ao(cost);    % thanks to markus buehren    
    c = 1;   % the index of the objects to be carried over
    clear objcarry;    % clear up carry over matrix
    
    % now assignment the specific objects to the tracks
    for m = 1:M       % loop through all the parental objects to be tracked
        n = assign(m);
        if (n)   % an assignment was made
            t1 = objnow(m).t;
            t2 = (i+1);
            ind1 = objnow(m).ind;
            ind2 = n;
            trno = objnow(m).tr;
            objects(t1).obj(ind1).trno = trno;   % assign the track number to the object in the first frame
            objects(t2).obj(ind2).trno = trno;   % assign the same track number to the object in second frame
            
            Ls(trno) = Ls(trno) + 1;   % increment to the length of the track
            
        else     % an assignment was NOT made, check to see if object needs to be carried over
            objnow(m).skipped = objnow(m).skipped + 1;   % increment the skipped frame number for this particular object
            if (objnow(m).skipped <= maxskip)     % carry this object over to the last frame
                objcarry(c) = objnow(m);
                c = c+1;
            end
        end
    end
end
save(out, 'objects', 'gatenames');

end
