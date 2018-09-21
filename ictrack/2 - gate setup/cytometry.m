function varargout = cytometry(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cytometry_OpeningFcn, ...
                   'gui_OutputFcn',  @cytometry_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

function cytometry_OpeningFcn(hObject, eventdata, handles, varargin)

%% load the image acquisition file
% mainly for getting the time information
[filename, dirname] = uigetfile('*.mat','Choose Acquisition file...');
load([dirname filename]);
handles.acq = acq;

%% find out which frames are common in all frames
segchannel = 2;  % segmentation channel, mCherry in this case
tseg = acq.C(segchannel).tlist;
tcommon = find(tseg);   % the list of all planes where all channels are present
thr = (acq.tr - acq.tr(1)).*24;   % time in hours
thr = thr(tcommon);     % this is the time in hours

%% load the cell segmentation file
[filename, dirname] = uigetfile('*.mat','Choose Schnitz file...');
load([dirname filename]);
T = length(objects);
withtrack = exist('tracks');

if withtrack  % initialize the Schnitz object using data from the file 
    sch = Schnitz2(objects, tracks);
else
    sch = Schnitz2(objects);
end
handles.sch = sch;

%% initialize data structure, X/Y popup menus
data = [];
for t = 1:T        
    obj = objects(t).obj;    
    data1 = cat(2,obj.data);
    delobj = [];
    for i = 1:length(data1)
        
        data1(i).t = thr(t);   % this is the time in hours
        data1(i).ind = i;
        data1(i).lgfp = log10(data1(i).gfp);
        data1(i).lrfp = log10(data1(i).rfp);                
        if (~withtrack)
            continue
        end        
        %% enter track information if present
        trno = obj(i).trno;
        
        if (trno == -1)
            delobj = [delobj i];
            continue;
        end
           
        if (trno)    
            data1(i).trno = trno;
            data1(i).approved = tracks(trno).approved;
        else
            data1(1).trno = 0;
            data1(i).approved = 0;
        end
    end    
    data1(delobj) = [];   % delete the deleted objects
    data = [data data1];    
end
names = fieldnames(data);
set(handles.X_popup,'String',names);
set(handles.Y_popup,'String',names);

%% compute data limits for each field
bins = 100;
for i = 1:length(names)
    fieldinfo(i).name = names{i};   % the name    
    data1 = [data.(names{i})];   % get the list of all
    dmin = min(data1);
    dmax = max(data1);    
    fieldinfo(i).range = dmin:(dmax-dmin)./bins:dmax;     % the display range
    
    if strmatch(names{i},'t')
        fprintf('Set time vector...\n');
        fieldinfo(i).range = thr;
    end    
end

%% initialize the root gate
gates(1).name = 'root';
gates(1).ingate = ones(1,length(data));  % every cell is in the root gate
gates(1).parent = 0;
gates(1).children = 0;

%% save data 
handles.data = data;
handles.fieldinfo = fieldinfo;
handles.objects = objects;
handles.gates = gates;
handles.output = hObject;
guidata(hObject, handles);

function varargout = cytometry_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function refresh_Callback(hObject, eventdata, handles)

data = handles.data;
fieldinfo = handles.fieldinfo;
gates = handles.gates;

Xf = get(handles.X_popup,'Value');    % get the field number for the X axis plot
Yf = get(handles.Y_popup,'Value');    % get the field numebr for the Y axis plot
G = get(handles.gates_select,'Value'); % get the currently selected gate number

X = fieldinfo(Xf).range;    % the limits for the field displayed on the X axis
Y = fieldinfo(Yf).range;    % the limits for the field displayed on the Y axis
inds = find(gates(G).ingate);  % indices of cells that are in the gate

%% now plot the 2D histogram
figure(10); hold off;
Xdata = [data.(fieldinfo(Xf).name)];
Ydata = [data.(fieldinfo(Yf).name)];
Xdata = Xdata(inds);
Ydata = Ydata(inds);
Z = hist2d([Ydata', Xdata'], Y, X);
imshow(Z,[],'XData',X,'YData',Y);
axis xy; axis on; axis normal; colormap jet;
xlabel(fieldinfo(Xf).name);
ylabel(fieldinfo(Yf).name);

%% plot the gate if one is selected
g = get(handles.gates_view,'Val');
if (g>1)    
    Xfg = gates(g).Xf;
    Yfg = gates(g).Yf;    
    if (Xfg == Xf)&(Yfg == Yf)   % correct plots are displayed
        hold on;
        plot(gates(g).xi, gates(g).yi,'Color','w','LineWidth', 1.5);
    end
end
   

%% save the new ranges
fieldinfo(Xf).range = X;
fieldinfo(Yf).range = Y;
handles.fieldinfo = fieldinfo;
guidata(hObject,handles);


function X_popup_Callback(hObject, eventdata, handles)

function X_popup_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Y_popup_Callback(hObject, eventdata, handles)

function Y_popup_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function X_limits_Callback(hObject, eventdata, handles)
fieldinfo = handles.fieldinfo;
Xf = get(handles.X_popup,'Val');  % the field number for X
fieldinfo(Xf).range = getrange(Xf, fieldinfo);
handles.fieldinfo = fieldinfo;
guidata(hObject, handles);

function Y_limits_Callback(hObject, eventdata, handles)
fieldinfo = handles.fieldinfo;
Xf = get(handles.Y_popup,'Val');  % the field number for X
fieldinfo(Xf).range = getrange(Xf, fieldinfo);
handles.fieldinfo = fieldinfo;
guidata(hObject, handles);

function range = getrange(f, fieldinfo)    
    range = fieldinfo(f).range;
    dmin = min(range);
    dmax = max(range);
    bins = length(range)-1;
    r = inputdlg(['Enter min (' num2str(dmin) '), max (' num2str(dmax) '), and bins (' num2str(bins)  ') for ' fieldinfo(f).name ':']);        
    r = r{1};    
    r = eval(r);    
    range = r(1):(r(2)-r(1))./r(3):r(2);
         
function new_rect_Callback(hObject, eventdata, handles)

%% create a new roi rectangle
rlims = inputdlg('Enter limits of rectangle: [x1 y1 W H]');
figure(10);
if isempty(rlims{1});    % draw rectangle interactively
    h = imrect;
    lims = h.getPosition;    
    
else
    lims = eval(rlims{1});    % draw rectangle using specified coordinates    
    h = imrect(gca, lims);
end
xi = [lims(1); lims(1); lims(1)+lims(3); lims(1)+lims(3); lims(1)];
yi = [lims(2); lims(2)+lims(4); lims(2)+lims(4); lims(2); lims(2)];

%% redraw the gate;
bw = h.createMask;
handles = new_gate(handles, bw, xi, yi);
hold on;
delete(h)
plot(xi,yi,'Color','w','LineWidth',2);
guidata(hObject,handles);


function handles = new_gate(handles, bw, xi, yi);

%% start to assign cells into gates
Xf = get(handles.X_popup,'Val');   % these are the identities of the X and Y axes
Yf = get(handles.Y_popup,'Val');   % these are the identities of the X and Y axes

%% load relevant information
gates = handles.gates;
fieldinfo = handles.fieldinfo;
data = handles.data;
parent = get(handles.gates_select,'Val');  % this is the currently used gate

%% return the X and Y bin numbers for all objects 
Xdata = [data.(fieldinfo(Xf).name)];
Ydata = [data.(fieldinfo(Yf).name)];
xmin = fieldinfo(Xf).range(1);
xmax = fieldinfo(Xf).range(end);
xB = length(fieldinfo(Xf).range);
ymin = fieldinfo(Yf).range(1);
ymax = fieldinfo(Yf).range(end);
yB = length(fieldinfo(Yf).range);
Xbnum = floor((xB-1).*(Xdata - xmin)./(xmax-xmin))+1;   % bin numbers for each data point
Ybnum = floor((yB-1).*(Ydata - ymin)./(ymax-ymin))+1; 
Xbnum(find(Xbnum<1)) = 1;
Xbnum(find(Xbnum>(xB-1))) = xB-1;
Ybnum(find(Ybnum<1)) = 1;
Ybnum(find(Ybnum>(yB-1))) = yB-1;

%% find cells within the gate
mask = find(bw);
linYX = sub2ind([yB-1,xB-1],Ybnum,Xbnum);
ingate = ismember(linYX,mask);
ingate = gates(parent).ingate.*ingate;   % this is the indices of all objects that fall inside this gate and also only in the parent gate 

%% make the name of the gate
gatename = inputdlg('Enter name of the gate:');
gatename = [gates(parent).name ':' gatename{1}];  % gatename, with parent identity appended

% display gate statistics
title([num2str(sum(ingate(:))) ' of ' num2str(length(ingate)) ' cells are in the gate.']);

%% save gate information
G = length(gates)+1;
gates(G).name = gatename;
gates(G).ingate = ingate;
gates(G).parent = parent;
gates(parent).children = 1;    % the parent gate now has a child
gates(G).children = 0;
gates(G).Xf = Xf;
gates(G).Yf = Yf;
gates(G).xi = xi;
gates(G).yi = yi;
gatenames = {gates.name};
set(handles.gates_select,'String',gatenames);
set(handles.gates_view,'String',gatenames);
handles.gates = gates;


function gates_select_Callback(hObject, eventdata, handles)

function gates_select_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function gates_view_Callback(hObject, eventdata, handles)

function gates_view_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function save_Callback(hObject, eventdata, handles)
[filename, dirname] = uiputfile('*.mat','Save fieldinfo and gates file...');
gates = handles.gates;
fieldinfo = handles.fieldinfo;
save([dirname '/' filename], 'gates','fieldinfo');

function load_Callback(hObject, eventdata, handles)
[filename, dirname] = uigetfile('*.mat','Load fieldinfo and gates file...');
load([dirname '/' filename]);
gatenames = {gates.name};
set(handles.gates_select,'String',gatenames);
set(handles.gates_view,'String',gatenames);
handles.gates = gates;
handles.fieldinfo = fieldinfo;
guidata(hObject,handles);

function partition_Callback(hObject, eventdata, handles)

gates = handles.gates;
data = handles.data;
objects = handles.objects;



%% obtain list of terminal gates
termgates = [gates.children];  % list of terminal gates
termgates = find(~termgates);  

%% verify that terminal gates are disjoint
% do this later

%% this is the vector of gate ids 
ids = zeros(1, length(data));
for i = 1:length(termgates)    
    ingate = gates(termgates(i)).ingate;    
    if (sum(ingate.*ids)>0)
        error('Gates are not disjoint!');
    end   
    ids = ids + i.*ingate;    
end
gatenames = {gates.name};
for i = 1:length(data)            
    t = data(i).t;
    ind = data(i).ind;    
    objects(t).obj(ind).gate = ids(i);
end
[filename, dirname] = uiputfile('*.mat','Enter name of Schnitz file...');
save([dirname '/' filename], 'objects','gatenames');

%function mHist = hist2d ([vY, vX], vYEdge, vXEdge)
%2 Dimensional Histogram
%Counts number of points in the bins defined by vYEdge, vXEdge.
%size(vX) == size(vY) == [n,1]
%size(mHist) == [length(vYEdge) -1, length(vXEdge) -1]
%
%EXAMPLE
%   mYX = rand(100,2);
%   vXEdge = linspace(0,1,10);
%   vYEdge = linspace(0,1,20);
%   mHist2d = hist2d(mYX,vYEdge,vXEdge);
%
%   nXBins = length(vXEdge);
%   nYBins = length(vYEdge);
%   vXLabel = 0.5*(vXEdge(1:(nXBins-1))+vXEdge(2:nXBins));
%   vYLabel = 0.5*(vYEdge(1:(nYBins-1))+vYEdge(2:nYBins));
%   pcolor(vXLabel, vYLabel,mHist2d); colorbar
function mHist = hist2d (mX, vYEdge, vXEdge)
nCol = size(mX, 2);
if nCol < 2
    error ('mX has less than two columns')
end

nRow = length (vYEdge)-1;
nCol = length (vXEdge)-1;

vRow = mX(:,1);
vCol = mX(:,2);

mHist = zeros(nRow,nCol);

for iRow = 1:nRow
    rRowLB = vYEdge(iRow);
    rRowUB = vYEdge(iRow+1);
    
    vColFound = vCol((vRow > rRowLB) & (vRow <= rRowUB));
    
    if (~isempty(vColFound))
        
        
        vFound = histc (vColFound, vXEdge);
        
        nFound = (length(vFound)-1);
        
        if (nFound ~= nCol)
            disp([nFound nCol])
            error ('hist2d error: Size Error')
        end
        
        [nRowFound, nColFound] = size (vFound);
        
        nRowFound = nRowFound - 1;
        nColFound = nColFound - 1;
        
        if nRowFound == nCol
            mHist(iRow, :)= vFound(1:nFound)';
        elseif nColFound == nCol
            mHist(iRow, :)= vFound(1:nFound);
        else
            error ('hist2d error: Size Error')
        end
    end    
end

function maketrack_Callback(hObject, eventdata, handles)

%% load schnitz data
sch = handles.sch;
gates = handles.gates;
data = handles.data;
acq = handles.acq;

segchannel = 2;  % segmentation channel, mCherry in this case
tcommon = find(acq.C(segchannel).tlist);   % the list of all planes where all channels are present
thr = handles.acq.tr(tcommon);

%% error checking
% -1 check for self-referntial parent child relationships
% can remove this later in the next version of the tracking
% for i = 1:length(sch.Tr)
%     children = sch.tracks(i).children;
%     if sum(ismember(i,children));
%         sch.tracks(i).children = [];
%         fprintf(['Track ' num2str(i) ' found to be self-referential.  Removing children.\n']);
%     end
% end
% sch = sch.checklins();
% sch = sch.approveall();   

%% process all approved tracks

g = get(handles.gates_select,'Val');
fprintf(['Processing tracks in gate: ' gates(g).name '\n']);
datag = data(find(gates(g).ingate));   % find all the data points within this one particular gate
tr_ingate = unique([datag.trno]);
ingate_count = [datag.trno];        % allows counting of how many data points within a given track fall within the gate 

%% determine fraction of object in parent gate
g2 = gates(g).parent;
datag2 = data(find(gates(g2).ingate));   % all the data points
ingate2_count = [datag2.trno];

%% determine which tracks are more than a certain fraction in the gate
ingate_real = zeros(size(tr_ingate));

for i = 1:length(tr_ingate)    
    overlap = length(find(ingate_count == tr_ingate(i)));
    total = length(find(ingate2_count == tr_ingate(i)));    
    ingate_real(i) = overlap/total;
end
alltracks = tr_ingate(find(ingate_real > 0.5));
alltracks = intersect(find([sch.approved]), alltracks);  % can modify this to incorporate tracks that simply fall into a particular gate
L = length(alltracks);   % the initial number of tracks to be processed 
i = 1;   % generate all the tracks

keyboard

while (~isempty(alltracks))
    trs = sch.getlintr(alltracks(1));  % get the tracks for each lineage        
    for j = 1:length(trs) % determine whether the trajectory falls within gate or not
        overlap = length(find(ingate_count == trs(j).tr));        
        if intersect(alltracks, trs(j).tr)
            trs(j).ingate = 1;
        else
            trs(j).ingate = 0;
        end        
%         total =  length(trs(j).ts);
%         if (overlap/total > 0.5)   % we just need trajectory to overlap more than 50% in gate
%             fprintf(['Track ' num2str(trs(j).tr) ' overlaps ' num2str(overlap/total) ' with gate.  Added\n']);
%             trs(j).ingate = 1;
%         else
%             fprintf(['Track ' num2str(trs(j).tr) ' overlaps ' num2str(overlap/total) ' with gate.  Not Added\n']);
%             trs(j).ingate = 0;
%         end                    
    end    
    [tracked, inda, indb] = intersect(alltracks, [trs.tr]);    
    alltracks(inda) = [];   % delete from the list all tracks that have included in these tracks    
    if (isempty(intersect(inda, 1)))
        fprintf('cannot retrace graph\n');
        pause;
    end          
    lin(i).trs = trs;         
    i = i+1;
    
end
keyboard
handles.tracked = 1;

%% save the lineage.mat file
[filename dirname] = uiputfile('lin.mat','Save lineages file');
save([dirname '/' filename], 'lin','thr');
handles.sch = sch;
guidata(hObject, handles);


function new_poly_Callback(hObject, eventdata, handles)
% create a new roipoly 
figure(10);
title('Choose gate...');
h = impoly;
pos = h.getPosition;
xi = [pos(:,1) ; pos(1,1)];
yi = [pos(:,2) ; pos(1,2)];

bw = h.createMask;
handles = new_gate(handles, bw, xi, yi);
hold on;
delete(h)
plot(xi,yi,'Color','w','LineWidth',2);
guidata(hObject,handles);
