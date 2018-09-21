function varargout = cytometry2(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cytometry2_OpeningFcn, ...
                   'gui_OutputFcn',  @cytometry2_OutputFcn, ...
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

function cytometry2_OpeningFcn(hObject, eventdata, handles, varargin)

%% import the image acquisition file
% mainly for getting the time information
[filename, dirname] = uigetfile('*.mat','Choose Acquisition file...');
load([dirname filename]);
handles.acq = acq;

%% find out which frames are common in all frames
thr = gettime(acq, 2);    % get times in hours where channel 2 (mCherry) is present

%% import the cell segmentation file, determine whether variables exist
[filename, dirname] = uigetfile('*.mat','Choose Schnitz file...');
load([dirname filename]);
withtrack = exist('tracks');
withgates = (exist('gates') && exist('fieldinfo'));

%% initialize Schnitz and Fcs classes, depending on what variables exist
if (withtrack) & (withgates)     % both gate and tracks exist                    
    sch = Schnitz2(timepoints, tracks);           
    fcs = Fcs(timepoints, thr, tracks, fieldinfo, gates)
    
elseif (withtrack) & (~withgates)    % track only
    sch = Schnitz2(timepoints, tracks);   
    fcs = Fcs(timepoints, thr, tracks)
    
elseif (~withtrack) & (withgates)    % gates only
    sch = Schnitz2(timepoints);
    fcs = Fcs(timepoints, thr, [], fieldinfo, gates)
    
else                           % only the raw objects exist
    sch = Schnitz2(timepoints); 
    fcs = Fcs(timepoints, thr)
end

%% initialize popup menus
set(handles.X_popup,'String',fcs.fnames);
set(handles.Y_popup,'String',fcs.fnames);
set(handles.gates_select,'String',fcs.gnames);
set(handles.gates_view,'String',fcs.gnames);

%% export the Schnitz and Fcs classes
handles.sch = sch;
handles.fcs = fcs;
handles.thr = thr;
handles.output = hObject;
guidata(hObject, handles);

function varargout = cytometry2_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function refresh_Callback(hObject, eventdata, handles)
%% import variables
fcs = handles.fcs;
data = handles.fcs.data;
fieldinfo = handles.fcs.fieldinfo;
gates = handles.fcs.gates;

%% get fields for display from the GUI
Xf = get(handles.X_popup,'Value');    % get the field number for the X axis plot
Yf = get(handles.Y_popup,'Value');    % get the field number for the Y axis plot
G = get(handles.gates_select,'Value'); % get the currently selected gate number
gv = get(handles.gates_view, 'Val');

%% retrieve data and data range
X = fcs.getrange(Xf);
Y = fcs.getrange(Yf);
Xdata = fcs.getdata(Xf, G);
Ydata = fcs.getdata(Yf, G);

Z = hist2d([Ydata, Xdata], Y, X);

%% now plot the 2D histogram
figure(10); hold off;

imshow(Z,[],'XData',X,'YData',Y);
axis xy; axis on; axis normal; 
colormap jet;
%xlabel(fcs.fnames{Xf});
%ylabel(fcs.fnames{Yf});

%% plot the 1D histogram 
figure(11); hold off;
hist(Xdata, X);



%% plot the gate if one is selected
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

Xf = get(handles.X_popup,'Val');  % the field number for X
handles.fcs = handles.fcs.setrange(Xf);
guidata(hObject, handles);

function Y_limits_Callback(hObject, eventdata, handles)
Yf = get(handles.Y_popup,'Val');  % the field number for X
handles.fcs = handles.fcs.setrange(Yf);
guidata(hObject, handles);

        
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

%% donwload class
fcs = handles.fcs;

%% start to assign cells into gates
Xf = get(handles.X_popup,'Val');   % these are the identities of the X and Y axes
Yf = get(handles.Y_popup,'Val');   % these are the identities of the X and Y axes
parent = get(handles.gates_select,'Val');  % this is the currently used gate

%% make the name of the gate
gatename = inputdlg('Enter name of the gate:');

fcs = fcs.newgate(gatename{1}, bw, xi, yi, Xf, Yf, parent);

set(handles.gates_select,'String',fcs.gnames);
set(handles.gates_view,'String',fcs.gnames);

handles.fcs = fcs;  % need to change this code for later


function gates_select_Callback(hObject, eventdata, handles)

function gates_select_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function gates_view_Callback(hObject, eventdata, handles)

gv = get(handles.gates_view, 'Val');
Xf = get(handles.X_popup,'Value');    % get the field number for the X axis plot
Yf = get(handles.Y_popup,'Value');    % get the field number for the Y axis plot
fcs = handles.fcs;

if (fcs.isgate(gv, Xf, Yf))
    hold on;    
    xi = fcs.gates(gv).xi;
    yi = fcs.gates(gv).yi;
    if ~isempty(xi)
        figure(10); hold on; plot(xi, yi,'Color','w','LineWidth', 1.5);
        figure(12); hold on; plot(xi, yi,'Color','k','LineWidth', 2);
    end
end

function gates_view_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function export_Callback(hObject, eventdata, handles)
[filename, dirname] = uiputfile('*.mat','Save fieldinfo and gates file...');
gates = handles.gates;
fieldinfo = handles.fieldinfo;
save([dirname '/' filename], 'gates','fieldinfo');

function import_Callback(hObject, eventdata, handles)
[filename, dirname] = uigetfile('*.mat','Load fieldinfo and gates file...');
load([dirname '/' filename]);
gatenames = {gates.name};
set(handles.gates_select,'String',gatenames);
set(handles.gates_view,'String',gatenames);
handles.gates = gates;
handles.fieldinfo = fieldinfo;
guidata(hObject,handles);

function save_Callback(hObject, eventdata, handles)

fieldinfo = handles.fcs.fieldinfo;
gates = handles.fcs.gates;
data = handles.fcs.data;
timepoints = handles.sch.timepoints;
tracks = handles.sch.tracks;

%% obtain list of terminal gates
termgates = [gates.children];  % list of terminal gates
termgates = find(~termgates);  

%% verify that terminal gates are disjoint
% do this later

%% this is the vector of gate ids 
ids = zeros(length(data), 1);
for i = 1:length(termgates)    
    ingate = gates(termgates(i)).ingate;    
    if (sum(ingate.*ids)>0)
        error('Gates are not disjoint!');
    end   
    ids = ids + i.*ingate;    
end
gatenames = {gates.name};

for i = 1:length(data)        
    t = data(i).tint;    % this is the time in frames
    ind = data(i).ind;
    timepoints(t).obj(ind).gate = ids(i);
end
[filename, dirname] = uiputfile('*.mat','Enter name of Schnitz file...');
save([dirname '/' filename], 'timepoints','tracks','gatenames','fieldinfo','gates');


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

%% import schnitz data
sch = handles.sch;
gates = handles.fcs.gates;
data = handles.fcs.data;
acq = handles.acq;

segchannel = 2;  % segmentation channel, mCherry in this case
tcommon = find(acq.C(segchannel).tlist);   % the list of all planes where all channels are present
thr = gettime(acq,segchannel);

%handles.acq.tr(tcommon);

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
handles.tracked = 1;

%% export the lineage.mat file
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


% --- Executes on button press in showmultigate.
function showmultigate_Callback(hObject, eventdata, handles)
%% import variables
fcs = handles.fcs;
data = handles.fcs.data;
fieldinfo = handles.fcs.fieldinfo;
gates = handles.fcs.gates;

pops = inputdlg('Enter the populations desired...');
pops = eval(pops{1});

Xf = get(handles.X_popup,'Value');    % get the field number for the X axis plot
Yf = get(handles.Y_popup,'Value');    % get the field number for the Y axis plot
%gv = get(handles.gates_view, 'Val');
X = fcs.getrange(Xf);
Y = fcs.getrange(Yf);

for i = 1:length(pops)
    %% retrieve data and data range
    G = pops(i);
    Xdata = fcs.getdata(Xf, G);
    Ydata = fcs.getdata(Yf, G);
    Z = hist2d([Ydata, Xdata], Y, X);
    d(i).Z = Z;
end

plot_multicolor(d,X,Y);

% keyboard
% 
% %% now plot the 2D histogram
% figure(10); hold off;
% imshow(Z,[],'XData',X,'YData',Y);
% axis xy; axis on; axis normal; 
% colormap jet;
% %xlabel(fcs.fnames{Xf});
% %ylabel(fcs.fnames{Yf});
% 
% %% plot the gate if one is selected
% if (fcs.isgate(gv, Xf, Yf))
%     hold on;
%     plot(fcs.gates(gv).xi, fcs.gates(gv).yi,'Color','w','LineWidth', 1.5);
% end

guidata(hObject,handles);

function loadlin_Callback(hObject, eventdata, handles)

fcs = handles.fcs;   % fcs object
Trf = strmatch('trno',fcs.fnames);   % this is the field number of the trno field
trnos = fcs.getdata(Trf);

[file dir] = uigetfile('*.mat', 'Enter name of lineage file...');
load([dir file]);

linno = zeros(size(trnos));
inlin = zeros(size(trnos));

% find data points corresponding to lineages and tracks in lineage file
for i = 1:length(lin)   % loop through all the lineages
    trs = lin(i).trs;

    for j = 1:length(trs)   % loop through all tracks in individual lineages        
        if isempty(trs(j).tr)
            continue   % if track has been deleted from lineage, then don't even list it
        end        
        trno = trs(j).tr;
        ingate = trs(j).ingate;
        inds = find(trnos == trno);   % find where all track data points are        
        linno(inds) = i;
        inlin(inds) = ingate;
    end
end

linn = inputdlg('Enter name of tracks...');
linn = linn{1};

fcs = fcs.newfield([linn '_linno'],linno);
fcs = fcs.newfield([linn '_inlin'],inlin);

set(handles.X_popup,'String',fcs.fnames);
set(handles.Y_popup,'String',fcs.fnames);

handles.fcs = fcs;
guidata(hObject, handles);

% need to add an additional field to the 