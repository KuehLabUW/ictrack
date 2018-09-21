function varargout = pview(varargin)
% Last Modified by GUIDE v2.5 13-Mar-2012 16:14:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @pview_OpeningFcn, ...
    'gui_OutputFcn',  @pview_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
end

function pview_OpeningFcn(hObject, eventdata, handles, varargin)

%% load acquisition file
[imfile imdir] = uigetfile('*.mat', 'Point to acq.mat file...');
load([imdir imfile]);   
if ~exist('acq')
    error('Cannot find ACQ acquisition structure in file...');
end
%acq.Z = 1; %length(acq.Z.bottom:acq.Z.step:acq.Z.top);   % update this in a future version of the pview
acq.tr = (acq.tr - acq.tr(1)).*24;   % convert from days to hours
acq.imdir = imdir;

%% load schnitz data
cd(imdir);
[imfile2 imdir] = uigetfile('*.mat',' Point to schnitz.mat file...');
load([imdir imfile2]);   


if (~exist('timepoints'))
    if exist('objects');
        %Ensuring compatibility with image structures named OBJECTS' 
        timepoints = objects ;
    else
        error('invalid schnitz file...');
    end
end

if exist('tracks')  % initialize the Schnitz object using data from the file 
    sch = Schnitz3(timepoints, tracks);
else
    sch = Schnitz3(timepoints);
end
handles.sch = sch;

%% if track info exists, append to 

%% if gatenames exist, then append them to the gates_select popup menu
if exist('gatenames')
    fprintf('found gatenames\n');   
else    
    gatenames = {'root:'};
end
set(handles.gates_select,'String',gatenames);


%% find out which frames are common in all frames
segchannel = 2;  % segmentation channel, mCherry in this case
tseg = acq.C(segchannel).tlist;
tcommon = find(tseg);   % the list of all planes where all channels are present

% compile list of previous frames for mismatch overlaying
tprev(1) = 1;
for i = 2:length(tseg)
    tprev(i) = tseg(i).*i + (1-tseg(i)).*tprev(i-1);
end
tprev = tprev';
acq.tcommon = tcommon;
acq.tprev = tprev;
Tc = length(tcommon);     % this is the number of planes where segmentation channel is present 
handles.acq = acq;

%% set the t and z limits for the sliders
set(handles.plane_slider,'Max',acq.T,'Sliderstep', [1/(acq.T-1),10/(acq.T-1)]);  % set time slider length
set(handles.planecommon,'Max',Tc,'Sliderstep', [1/(Tc-1),10/(Tc-1)]);  % set time slider length


%keyboard
%if (acq.Z > 1)
%    set(handles.zslider,'Max',acq.Z,'Sliderstep', [1/(acq.Z-1),10/(acq.Z-1)]);  % set z slider length
%end

%% set the channel displays and the contrast values, if present
for i = 1:length(acq.C);
    set(handles.(['show' num2str(i)]), 'Value', 1);  % turn on the slider 
    if isfield(acq.C, 'min');
        set(handles.(['clow' num2str(i)]),'String',num2str(acq.C(i).min));
    end
    
    if isfield(acq.C, 'min');
        set(handles.(['chi' num2str(i)]),'String',num2str(acq.C(i).max));
    end
end

%% display first frame
handles = initfig(handles,1);  % initialize figure axes
handles = refreshfig(handles,1);     % refresh the figure, update the image buffer if necessary

%% save variables to handles structure
handles.tracked = 0;
guidata(hObject, handles);
end

% --- Outputs from this function are returned to the command line.
function varargout = pview_OutputFcn(hObject, eventdata, handles)
end

%%%%% MY OWN FUNCTIONS %%%%%
function handles = initfig(handles,fignum)   % initialize the figure display for a given offset
acq = handles.acq;
X = acq.X;
Y = acq.Y;
N = acq.N;
M = acq.M;
Xo = str2num(get(handles.xoffset,'String'));
Yo = str2num(get(handles.yoffset,'String'));

clf;   % initialize figure;

figure(1);
axis image; axis ij; axis off;   % initialize the figure axis
set(gcf, 'WindowButtonMotionFcn',{@updatetitle,handles});
set(gcf, 'MenuBar','none');
set(gcf, 'ToolBar','none');

axis([1 N*(X-Xo)+Xo 1 M*(Y-Yo)+Yo]);
set(gca,'Position', [0 0 1 1]);
colormap gray
hold on;

[meshx, meshy] = meshgrid(1:X, 1:Y);   % image indices
inx = NaN.*zeros(M*(Y-Yo)+Yo, N*(X-Xo)+Xo,M,N);
iny = inx;
% calculate the closest frame from a given point
closest = zeros(M*(Y-Yo)+Yo, N*(X-Xo)+Xo, M, N);
closem = zeros(size(closest,1),1);   % vector giving the closest m mindex
closen = zeros(size(closest,2),1);

for m = 1:M   % map closest row positions of images in the tile
    mmin = round(Yo/2 + (m-1)*(Y-Yo));
    mmax = round(Yo/2 + m*(Y-Yo));
    if (m==1)
        closem(1:mmax) = m;
    elseif (m == M)
        closem(mmin:end) = m;
    else
        closem(mmin:mmax) = m;
    end
end

for n = 1:N % map closest column positions of images in the tile
    nmin = round(Xo/2 + (n-1)*(X-Xo));
    nmax = round(Xo/2 + n*(X-Xo));
    
    if (n==1)
        closen(1:nmax) = n;
    elseif (n==N)
        closen(nmin:end) = n;
    else
        closen(nmin:nmax) = n;
    end
end

% map index in image tile to index in individual images, for segmentation purposes
% visible border

bvis = [];
binvis = [];

for m = 1:M
    for n = 1:N
        xmin = (n-1)*(X-Xo)+1;
        ymin = (m-1)*(Y-Yo)+1;
        xmax = xmin+X-1;
        ymax = ymin+Y-1;
        loc(m,n) = imagesc(xmin, ymin, zeros(Y, X));  % initialize handles for graphics display
        inx(ymin:ymax, xmin:xmax,m,n) = meshx;
        iny(ymin:ymax, xmin:xmax,m,n) = meshy;
        xmins(m,n) = xmin;
        xmaxs(m,n) = xmax;
        ymins(m,n) = ymin;
        ymaxs(m,n) = ymax;
        % the border of the image that is 'visible', i.e. top and left
        bvis = [bvis [xmin xmin xmax NaN ; ymax ymin ymin NaN]];   % xy positiony of upper leftborders
        binvis = [binvis [xmin xmax xmax NaN ; ymax ymax ymin NaN]];   % xy positions lower right borders
    end
end
tile.bvh = 0;   % handle for the visible border, zero if not shown
tile.bih = 0;
tile.segh = [];   % handle for the segmentations
tile.loc = loc;
tile.inx = inx;
tile.iny = iny;
tile.closem = closem;
tile.closen = closen;

tile.xmins = xmins;
tile.xmaxs = xmaxs;
tile.ymins = ymins;
tile.ymaxs = ymaxs;

tile.bvis = bvis;
tile.binvis = binvis;

% show the border if the showborder check box is checked
if (get(handles.showborder,'Value'))
    tile = makeborder(tile);  % function for displaying the border and storing the resultant handle in the tile structure;
end
handles.tile = tile;   % save all the tiling information back to the handles structure

% set the elapsed time indicator
handles.timeh = text(0.85, 0.06,'time','Units','Normalized','Color','w','FontUnits','Normalized','FontSize', 0.06);
end

function tile = makeborder(tile)
figure(1);
tile.bvh = plot(tile.bvis(1,:), tile.bvis(2,:), 'y', 'LineWidth',2);
%tile.bih = plot(tile.binvis(1,:), tile.binvis(2,:), 'y:', 'LineWidth',1);
end

function handles = refreshfig(handles,reload)     % refresh the figure, update the image buffer if necessary

tnow = str2num(get(handles.plane_current,'String'));  % the current frame
treal = handles.acq.tr(tnow);  % this is the real time
T = handles.acq.T;    % the total number of frames in the stack
M = handles.acq.M;
N = handles.acq.N;

%% display tiled image
if (reload)
    buffer(T).imall = [];  % clear the image buffer
    handles.buffer = buffer;
end
imall = handles.buffer(tnow).imall;

%fprintf(['time point ' num2str(tnow) '.image isempty flag is ' num2str(isempty(imall)) '...\n']);
if isempty(imall);
    imall = imload(handles);   % generate figure
    handles.buffer(tnow).imall = imall;   % store figure into buffer    
end

showfig(imall, handles.tile.loc);   % show figure

%% display elapsed time
if (get(handles.showtime,'Value'));    
    treal = round(treal*100)./100; % round to two decimal places
    
    set(handles.timeh,'String',[num2str(treal) '''']);
else
    set(handles.timeh,'String','');
end

%% display border of segmented and tracked objects
eraseline(handles.tile);   % erase all the line objects in figure 1    
tdiff = (tnow - handles.acq.tprev(tnow));  % this is the difference between the current frame and the last frame where segmentation was performed.//
if ((get(handles.showtrack,'Value'))&(tdiff==0))       
    % show all tracks in a single lineage
    sch = handles.sch;   % this is the schnitz information    
    tr_s = str2num(get(handles.tr,'String'));  % this is the selected track    
    tcnow = round(get(handles.planecommon,'Value'));
    obj = sch.timepoints(tcnow).obj;   % load the objects structure for the current time point
     
    if isempty(obj)
        return
    end    
    %% display only cells within a certain gate, if a gate is selected    
    
    cellgates = [obj.gate];    
    selectedgate = get(handles.gates_select,'Value') - 1;    % start numbering at zero

    
    
    %% display only objects within the x and y axis limits to save some time            
    
    tile = handles.tile;
    xminlist = tile.xmins(:);
    yminlist = tile.ymins(:);    
    locinds = sub2ind(size(tile.xmins),[obj.m]',[obj.n]');    
    xlist = xminlist(locinds) + [obj.x]';
    ylist = yminlist(locinds) + [obj.y]';    
    xlims = get(gca,'XLim');
    ylims = get(gca,'YLim');    
    % find only the objects that satisfy all the limits    
    onaxis = (xlist > xlims(1))&(xlist < xlims(2))&(ylist > ylims(1))&(ylist < ylims(2));
    onaxis = onaxis';
    
    if (~selectedgate)  % select all possible cells
        displayedcells = find(onaxis);        
    else
        displayedcells = find((cellgates == selectedgate)&onaxis);
    end
        
    for i = displayedcells   % cycle through all the objects        
        if (obj(i).trno > 0) % object is a tracked object                  
            rgb = sch.rgb(obj(i).trno);   % color settings the display            
            if (sch.approved(obj(i).trno))            
                style = ':';    % dotted line for object in approved track
            else
                style = '-';    % solid line for object in unapproved track
            end            
            if (obj(i).trno == tr_s)
                width = 2;
            else
                width = 1;
            end
            showobj(obj(i), handles.tile, rgb, width, style);
        elseif (obj(i).trno == 0);  % object is not tracked
            rgb = [1 1 1];
            style = '--';     % white dashed line for a unassigned object
            width = 1;
            showobj(obj(i), handles.tile, rgb, width, style);
        end        
    end
end
end

function [imall, varargout] = imload(handles)

showmissing = get(handles.showmissing,'Value');  % show missing frames in the image
acq = handles.acq;    % load the acquisition structure
tnow = str2num(get(handles.plane_current,'String'));  % load the current t and z slice
znow = str2num(get(handles.zcurrent,'String'));
% load contrast settings for individual channels
minc(1) = str2num(get(handles.clow1,'String'));
maxc(1) = str2num(get(handles.chi1,'String'));
minc(2) = str2num(get(handles.clow2,'String'));
maxc(2) = str2num(get(handles.chi2,'String'));
minc(3) = str2num(get(handles.clow3,'String'));
maxc(3) = str2num(get(handles.chi3,'String'));
minc(4) = str2num(get(handles.clow4,'String'));
maxc(4) = str2num(get(handles.chi4,'String'));
minc(5) = str2num(get(handles.clow5,'String'));
maxc(5) = str2num(get(handles.chi5,'String'));
minc(6) = str2num(get(handles.clow6,'String'));
maxc(6) = str2num(get(handles.chi6,'String'));
% load the show settings for individual channels
show(1) = get(handles.show1,'Value');
show(2) = get(handles.show2,'Value');
show(3) = get(handles.show3,'Value');
show(4) = get(handles.show4,'Value');
show(5) = get(handles.show5,'Value');
show(6) = get(handles.show6,'Value');

M = acq.M;   % acquisition file information
N = acq.N;
Z = acq.Z;
C = acq.C;
T = acq.T;
X = acq.X;
Y = acq.Y;
Z = 1; %acq.Z;
tprev = acq.tprev;  % this is the list of previous common images for showing missing frames
imdir = acq.imdir;
zproj = get(handles.zproj,'Value'); % do a z projection or not
images = loadone(imdir,tnow);   % load the image at the current frame
if (showmissing & (tprev(tnow) < tnow))   % if there are missing frames and previous images need to be loaded
    imprev = loadone(imdir,tprev(tnow));
end

imall = uint16(zeros(Y,X,3,M,N));   % first define the image i,j,color,M,N

for m = 1:M
    for n = 1:N
        imrgb = uint16(zeros(Y,X,3));
        for c = 1:length(C);
            existtnow = C(c).tlist(tnow);  % flag for whether channel exists at this particular time point
            if (~show(c))   % skip if not showing this channel altogether
                continue;   % go to he
            elseif ((~existtnow) & (~showmissing));  % images don't exist, and not showing missing frames
                continue;
            else % show an image
                if (existtnow)  % current time point exists
                    imuse = images;
                else
                    imuse = imprev;
                end
                % add the channel to the multi-color overlay
                if (C(c).doz)
                    if (zproj)    % sum of zs
                        im = imadjust16(imuse(c).im(:,:,m,n,:), Z*[minc(c) maxc(c)]);  % scaled image
                        im = uint16(sum(im,5));   % sum the image across all z slices
                    else
                        im = imadjust16(imuse(c).im(:,:,m,n,znow), [minc(c) maxc(c)]);  % scaled image
                    end
                else
                    im = imadjust16(imuse(c).im(:,:,m,n,1), [minc(c) maxc(c)]);  % scaled image
                end
                imrgb = imrgb + cat(3, im*(C(c).cR/255), im*(C(c).cG/255), im*(C(c).cB/255));
            end
        end
        %         if (length(images) >= 10)   % labels and segmentation exists
        %             labels(:,:,m,n) = images(10).im(:,:,m,n);   % load the labels for the current frame
        %         else                                    % load the labels for the previous frame
        %             labels(:,:,m,n) = [];
        %         end
        imall(:,:,:,m,n) = imrgb;
    end
end
% if (nargout == 2)
%    varargout(1) = {labels};
% end
end

function showfig(imall,loc)        % displays the image stk
for m = 1:size(loc,1)
    for n = 1:size(loc,2)
        set(loc(m,n),'CData', imall(:,:,:,m,n));
    end
end
drawnow;
end


function updatetitle(src,eventdata,handles);
currentxy = get(gca,'CurrentPoint');
xnow = round(currentxy(1,1));
ynow = round(currentxy(1,2));
int = 1;
set(src,'name', ['x = ' num2str(xnow) '; y = ' num2str(ynow) '; int = ' num2str(int)]);
end
%%%%%%%% Callback functions
function planecommon_Callback(hObject, eventdata, handles)
tcommon = handles.acq.tcommon;    % the list of all planes that have all channels
tcnow = round(get(hObject, 'Value'));   % get the current value of the slider
set(handles.plane_current,'String',num2str(tcommon(tcnow)));
set(handles.plane_slider,'Value', tcommon(tcnow));
handles = refreshfig(handles,0);
guidata(hObject,handles);
end

function planecommon_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

function plane_slider_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
set(hObject,'Value',1);
end

function plane_slider_Callback(hObject, eventdata, handles)
handles = setslider(handles);
handles = refreshfig(handles,0);
guidata(hObject, handles);
end

function handles = setslider(handles);
tnow = round(get(handles.plane_slider,'Value'));  % get the frame number of the current stack
set(handles.plane_current,'String', num2str(tnow));
tcommon = handles.acq.tcommon;  % list of all planes that have all channels
under = find(tnow-tcommon >= 0);  % find all frames with a frame number greater than
tcnow = under(end);
if isempty(tcnow)
    error('unknown error...');
else
    set(handles.planecommon, 'Value', tcnow);
end
end


function loadall_Callback(hObject, eventdata, handles)

T = handles.acq.T;   % all the planes in the movie
tnow = str2num(get(handles.plane_current,'String'));

for i = 1:T
    set(handles.plane_slider,'Value',i);
    set(handles.plane_current,'String',num2str(i));
    if (i==1)
        handles = refreshfig(handles,1);  % set the reload
    else
        handles = refreshfig(handles,0);  % do not reload,
    end
end
% set the plane back to the previous plane before the program started
set(handles.plane_slider,'Value',tnow);
set(handles.plane_current,'String',num2str(tnow));
handles = refreshfig(handles,0);
guidata(hObject, handles);
end

function saveimage_Callback(hObject, eventdata, handles)
% get the file
[filename, pathname] = uiputfile('*.tif', 'Save image file as...')

if isempty(filename)
    return
end

figure(1);
v = axis;
imin = ceil(v(3));
imax = floor(v(4));
jmin = ceil(v(1));
jmax = floor(v(2));

t = handles.t;
im = handles.ims(t).im;
subim = im(imin:imax, jmin:jmax, :);
imwrite(subim, [pathname, filename], 'Compression', 'None');
figure(handles.figure1);
end

%%%%%% accessory callback functions %%%%%%%
function cutregion_Callback(hObject, eventdata, handles)
end

function bridgeregion_Callback(hObject, eventdata, handles)
end

function show1_Callback(hObject, eventdata, handles)
handles = refreshfig(handles,1);
guidata(hObject, handles);
end

function clow1_Callback(hObject, eventdata, handles)
handles = refreshfig(handles,1);
guidata(hObject, handles);
end

function clow1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function chi1_Callback(hObject, eventdata, handles)
handles = refreshfig(handles,1);
guidata(hObject, handles);
end

function chi1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function show2_Callback(hObject, eventdata, handles)
handles = refreshfig(handles,1);
guidata(hObject, handles)
end

function clow2_Callback(hObject, eventdata, handles)
handles = refreshfig(handles,1);
guidata(hObject, handles);
end


function clow2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function chi2_Callback(hObject, eventdata, handles)
handles = refreshfig(handles,1);
guidata(hObject, handles);
end


function chi2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function show3_Callback(hObject, eventdata, handles)
handles = refreshfig(handles,1);
guidata(hObject, handles);
end

function clow3_Callback(hObject, eventdata, handles)
handles = refreshfig(handles,1);
guidata(hObject, handles);
end


function clow3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function chi3_Callback(hObject, eventdata, handles)
handles = refreshfig(handles,1);
guidata(hObject, handles);
end

function chi3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function show4_Callback(hObject, eventdata, handles)
handles = refreshfig(handles,1);
guidata(hObject, handles);
end


function clow4_Callback(hObject, eventdata, handles)
handles = refreshfig(handles,1);
guidata(hObject, handles);
end

function clow4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function chi4_Callback(hObject, eventdata, handles)
handles = refreshfig(handles,1);
guidata(hObject, handles);
end

function chi4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function show5_Callback(hObject, eventdata, handles)
handles = refreshfig(handles,1);
guidata(hObject, handles);
end
function clow5_Callback(hObject, eventdata, handles)
handles = refreshfig(handles,1);
guidata(hObject, handles);
end

function clow5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function chi5_Callback(hObject, eventdata, handles)
handles = refreshfig(handles,1);
guidata(hObject, handles);
end

function chi5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function show6_Callback(hObject, eventdata, handles)
handles = refreshfig(handles,1);
guidata(hObject, handles);
end

function clow6_Callback(hObject, eventdata, handles)
handles = refreshfig(handles,1);
guidata(hObject, handles);
end

function clow6_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function chi6_Callback(hObject, eventdata, handles)
handles = refreshfig(handles,1);
guidata(hObject, handles);
end

function chi6_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function zslider_Callback(hObject, eventdata, handles)
znow = round(get(hObject,'Value'));  % get the frame number of the current stack
set(handles.zcurrent,'String', num2str(znow));
handles = refreshfig(handles,1);
guidata(hObject, handles);
end

function zslider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

function zcurrent_Callback(hObject, eventdata, handles)
znow = str2num(get(hObject,'String'));
Z = handles.acq.Z;
if (znow > Z)  % error checking
    znow = Z;
elseif (znow < 1)
    znow = 1;
end
set(hObject,'String', num2str(znow));
set(handles.zslider,'Value',znow);
handles = refreshfig(handles,1);
guidata(hObject, handles);
end


function zcurrent_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function zproj_Callback(hObject, eventdata, handles)
handles = refreshfig(handles,1);
guidata(hObject,handles);
end
function xoffset_Callback(hObject, eventdata, handles)
handles = initfig(handles,1);
handles = refreshfig(handles,0);
guidata(hObject,handles);
end
function xoffset_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function yoffset_Callback(hObject, eventdata, handles)
handles = initfig(handles,1);
handles = refreshfig(handles,0);
guidata(hObject,handles);
end
function yoffset_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function trackno_Callback(hObject, eventdata, handles)

trackno = str2num(get(hObject,'String'));  % this is the current track that that is currently input
totaltracks = length(handles.tracks);

if (trackno > totaltracks)   % number exceeds the total track number
    fprintf('Track number exceeds the total number of tracks!\n');
    trackno = totaltracks;
end
% update the linno text box
linno = handles.tracks(trackno).lin;
set(handles.linno,'String',num2str(linno));
end

function trackno_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function showtrack_Callback(hObject, eventdata, handles)
handles = refreshfig(handles,0);
if (get(hObject,'Value')==0);
    eraseline(handles.tile);
end
guidata(hObject, handles);
end

function figure1_KeyPressFcn(hObject, eventdata, handles)
end

function segment_Callback(hObject, eventdata, handles)

% prevent stray keypress
set(gcf,'WindowStyle','modal');

set(gcf,'KeyPressFcn',@empty);
    function empty(src,eventdata);
    % empty function to prevent change of focus
    end

% basic image information
X = handles.acq.X;           
Y = handles.acq.Y;
T = handles.acq.T;
M = handles.acq.M;
maxzoom = X*M;  % assuming equal X and Y axis, if not, update code.

tcommon = handles.acq.tcommon;
Tc = length(tcommon);

% must turn on show tracks before any tracking can be done
set(handles.showtrack,'Value',1);  

% define keystrokes that require tracking and object selection
 
need_s = {'k','.','0'};   % needs selection but not track
need_tr = {'2','1', '3', '*', 'a', 's', 'x', '5', '\', '/', 'z','v','d','p'}; % needs track
need_tr_s = {'5', '\', '/', 'z'};   % needs track and selection

% Events structure to save actions
ev.save = 0;      % whether to save schnitz structure after a keypress
ev.tr = str2num(get(handles.tr,'String'));  % current selected track
ev.quit = 0;      % quit the 
ev.direction = 1; % 0 = no direction, 1 = right, -1 = left
ev.display = 'Tracking mode.  Press a button...';
ev.t = get(handles.plane_slider,'Value');
ev.tc = round(get(handles.planecommon,'Value'));
ev.undo = 0;

% copy of the schnitz structure.  save when needed
sch = handles.sch;    

while (~ev.quit)       
    
    % save the copy of sch into handles structure 
    if (ev.save)
        sch_old = handles.sch;   % for undo purposes     
        handles.sch = sch;
        ev.save = 0;
    end
    
    % update figure and text bar
    figure(1);  % refresh figure
    handles = refreshfig(handles,0); 
    set(handles.figure1,'Name',ev.display);   % display the text to be displayed
        
    % wait for keyboard button to be pressed 
    w = waitforbuttonpress; % get the next keypress (ignore if a mouse button was pressed)            
    inkey = get(gcf,'CurrentCharacter');
    
    switch inkey  % Master switch structure
            
        case '~'   % switch to the root gate 
            set(handles.gates_select,'Value',1);
            refreshfig(handles,0);
            
        case '!'   % switch to the first gate
            set(handles.gates_select,'Value',2);          
            refreshfig(handles,0);
            
        case '@'   % switch to the second gate
            set(handles.gates_select,'Value',3);            
            refreshfig(handles,0);
            
        case '#'   % switch to the third gate
            set(handles.gates_select,'Value',4);            
            refreshfig(handles,0);                   
    
        case '4' % go to the previous frame in the movie         
            ev.t = max(1, ev.t-1);
            set(handles.plane_slider,'Value',ev.t);            
            handles = setslider(handles);
            ev.tc = round(get(handles.planecommon,'Value'));   % get the current common timepoint                            

        case '6'  % go to the next frame in the movie            
            ev.t = min(ev.t+1,T);
            set(handles.plane_slider,'Value',ev.t);            
            handles = setslider(handles);
            ev.tc = round(get(handles.planecommon,'Value'));   % get the current common timepoint                
           

        case '7' % previous common frame                          
            ev.tc = max(1,round(get(handles.planecommon,'Value'))-1);
            ev.t = tcommon(ev.tc);
            set(handles.plane_slider,'Value',ev.t);
            handles = setslider(handles);            

        case '9' % next common frame    
            ev.tc = min(Tc, round(get(handles.planecommon,'Value'))+1);
            ev.t = tcommon(ev.tc);
            set(handles.plane_slider,'Value',ev.t);
            handles = setslider(handles);
            
        case '['  % set tracking direction to the left
            ev.direction = -1;
            ev.display = 'Tracking direction: backwards';
            
        case ']'  % set tracking direction to the right
            ev.direction = 1;
            ev.display = 'Tracking direction: forwards';
            
        case 'g' % go to given track
            tr = inputdlg('Enter track number to go to:');
            tr = str2num(tr{1});
            if (sch.istrack(tr))
                ev.tr = tr;
                tr_set(handles, ev.tr, sch.rgb(ev.tr));
                ev.display = ['Selected track ' num2str(ev.tr)];
            else
                ev.display = 'Track not found';
            end
            
        case '+'  % zoom
            myzoom(1.5, maxzoom);
            
        case '8' % pan
            myzoom(1, maxzoom);
            
        case '-' % zoom
            myzoom(1/1.5, maxzoom);            

        case '`'
            keyboard;
        
        case 'u' % undo
            if ~exist('sch_old')
                ev.display = 'Nothing to undo';
            else                
                sch = sch_old;               
                ev.save = 1;
                if (~ev.undo)
                    ev.display = 'Undo!';
                else
                    ev.display = 'Redo!';
                end
                ev.undo = ~ev.undo;
            end
            
        case 'q'  % quit            
            ev.quit = 1;            
            
        case need_s    % needs a selection
            [x, y] = gcaxy(gca);
            ind = identifyobj(handles, x, y);   % find index of the pointed object
            if (~ind)
                set(handles.figure1,'Name','No object found in frame');
                continue
            end
            
            switch inkey
                case '0'  % select current object, creating new track if necessary                    
                    if (~sch.tr(ev.tc,ind))   % object not already assigned to a track                        
                        sch = sch.addtotrack(0,ev.tc,ind);   % generate new track
                        ev.save = 1;
                    end                    
                    ev.tr = sch.tr(ev.tc,ind);                    
                    tr_set(handles, ev.tr, sch.rgb(ev.tr));                    
                    ev.display = ['Selected track ' num2str(ev.tr)];
                    
                case '.' % delete object                    
                    tr = sch.tr(ev.tc,ind);
                    if (~tr)   % object not assigned to track
                        sch = sch.delobj(ev.tc,ind);
                        ev.display = 'Object killed';
                        ev.save = 1;
                    else
                        % clear the selection if previous selected object is to be
                        % removed
                        ev.display = ['Object removed from track ' num2str(tr)];
                        sch = sch.removefromtrack(tr,ev.tc,ind);
                        ev.save = 1;
                    end
                    
                case 'k'   % kill entire track and all associated objects
                    % NEED TO MODIFY CODE TO KILL TRACK WITH MOUSE OVER
                    tr = sch.tr(ev.tc,ind);
                    if (~tr)   % object not assigned to track
                        sch = sch.delobj(ev.tc,ind);
                        ev.display = 'Object killed';
                        ev.save = 1;
                    elseif (tr == ev.tr)
                        ev.display = 'Cannot kill selected track.'
                    else                        
                        [ts inds] = sch.trackobj(tr);
                        for i = 1:length(ts)
                            sch = sch.delobj(ts(i), inds(i));
                        end
                        ev.save = 1;
                        ev.display = ['Track ' num2str(tr) ' killed.'];                                                
                    end
            end
            
%% these commands require a pre-existing track
        case need_tr   % the keystrokes require a track selection
            
            if (~ev.tr)
                ev.display = 'No track selected';
                continue
            end            
            
            switch inkey                        
                case '1' % go to the beginning of the current track                                                            
                    [ts, os] = sch.trackobj(ev.tr);
                    ev.tc = min(ts);
                    ev.t = tcommon(ev.tc);
                    set(handles.plane_slider,'Value',ev.t);
                    handles = setslider(handles);

                case '3' % go to the end of the current track                    
                    [ts, os] = sch.trackobj(ev.tr);
                    ev.tc = max(ts);
                    ev.t = tcommon(ev.tc);
                    set(handles.plane_slider,'Value',ev.t);
                    handles = setslider(handles);

               
                case '2'  % identify location of the current object                    
                    objs = sch.obj(ev.tc);
                    ind = find([objs.trno]==ev.tr);  % index of the object corresponding to the current track number                    
                    if (~isempty(ind))
                        m = objs(ind).m;
                        n = objs(ind).n;
                        x = objs(ind).x;
                        y = objs(ind).y;
                        %% code to display something                        
                        h = handles.tile.loc(m, n);  % these are the handles to the current image
                        xmin = get(h,'XData');  % access the x coordinate of the tiled sub-image
                        ymin = get(h,'YData');
                        plot(xmin+x-1, ymin+y-1,'Color','w','Marker','o','MarkerSize',12,'LineWidth',2);
                        ev.display = 'Identified object';
                    else
                        ev.display = 'Object not present in frame';                        
                    end
                    
%                 case 'v'  % view track lineage information
%                     set(handles.figure1,'Name','Making lineage tree....');                    
%                     thr = handles.acq.tr(tcommon(1:Tc));                    
%                     lintr = sch.getlintr(ev.tr);   % obtain the lineage information
%                     figure(3);
%                     subplot(1,3,1);
%                     showdata = get(handles.objdata,'Value');   % show the intensity data or not                    
%                     if showdata                        
%                         plotlin(thr, lintr, 'rfp','r',500000);
%                         xlabel('hours');
%                         set(gca,'XLim',[0 150]);
%                         set(gca,'YLim',[-1.2 1.2]);
%                         title('mCherry integrated intensity');
%                     else
%                         plotlin(thr, lintr);
%                     end
%                     
%                     subplot(1,3,2);                
%                     showdata = get(handles.objdata,'Value');   % show the intensity data or not                    
%                     if showdata                        
%                         plotlin(thr, lintr, 'gfp','g',200000);
%                         xlabel('hours');
%                         set(gca,'XLim',[0 150]);
%                         set(gca,'YLim',[-1.2 1.2]);
%                         title('PU.1-GFP integrated intensity');
%                     else
%                         plotlin(thr, lintr);
%                     end
%                     
%                     subplot(1,3,3);                
%                     showdata = get(handles.objdata,'Value');   % show the intensity data or not                    
%                     if showdata                        
%                         plotlin(thr, lintr, 'area',[0.5 0.5 0.5],1600);
%                         xlabel('hours');
%                         set(gca,'XLim',[0 150]);
%                         set(gca,'YLim',[-1.2 1.2]);
%                         title('Cell Area');
%                     else
%                         plotlin(thr, lintr);
%                     end
%                     
                    
                case 'p'   % draw a polygon to create a new object                                     
                    set(handles.figure1,'Name','Draw polygon (within single tile)');
                    figure(1);
                    V = [];
                    
                    while (isempty(V))|(V<3);   % perform if no polygon, or if number of vertices not enough 
                        h = impoly(gca);   % draw polygon
                        vs = round(getPosition(h));   % get the positions of all the vertices in the polygon
                        delete(h);        % delete polygon
                        V = size(vs,1);   % number of vertices
                    end
                    
                    v2s = [];                    
                    ms = [];
                    ns = [];
                    
                    % estimate centroid using mean x and y positions for
                    % entire tiled image
                    xmean = round(mean(vs(:,1)));
                    ymean = round(mean(vs(:,2))); 
                                        
                    for i = 1:size(vs,1)   % go through all the individual vertices in the polygon                        
                        [m,n,x,y] = tile2im_xy(handles.tile, vs(i,1), vs(i,2));                        
                        v2s = [v2s ; x y];                        
                        ms = [ms m];
                        ns = [ns n];
                    end                                            
                    
                    if (sum(abs(diff(ms))) | sum(abs(diff(ns))))    % vertices are not all on the same tile
                        error('vertices are not all on the same tile');
                    end                        
                    imb = poly2mask(v2s(:,1), v2s(:,2), Y, X);    % create mask of new object                    
                    objout = newobj(imb,m,n);   % generate new object                    
                    sch = sch.addobj(ev.tc,objout);   % add new objects
                    sch_old = handles.sch;        % manually save sch for undo
                    handles.sch = sch;              % need to update to allow new objects to be identified                            
                    ind2 = identifyobj(handles, xmean, ymean);     % obtain the index of the new object                    
                    sch = sch.addtotrack(ev.tr, ev.tc, ind2);
                    handles.sch = sch;
                    ev.save = 0;                  % performed this save manually to enable undo                            
                    ev.display = 'Generated new object';    
                    
                    % now need to make sure whether old object is handled
                    % correctly
                
                    
                    
                case '*'   % split track at the timepoint
                    sch = sch.splittrack(ev.tr, ev.tc+1);
                    ev.display = ['Split track ' num2str(ev.tr)];
                    ev.save = 1;
                    
                case 'a' % approve current track                    
                    sch = sch.approvetrack(ev.tr);
                    ev.display = ['Approved track # ' num2str(ev.tr)];
                    ev.save = 1;                    
                    
                case 's' % disapprove current track                     
                    sch = sch.approvetrack(ev.tr, 0);
                    ev.display = ['Disapproved track # ' num2str(ev.tr)];
                    ev.save = 1;
                    
                case 'x' % remove all children for given track                    
                    sch = sch.delchild(ev.tr);
                    ev.display = ['Removed children of track ' num2str(ev.tr)];
                    ev.save = 1;
                    
                case 'd'  % death flag on mark
                    data.death = 1;
                    sch = sch.addtrdata(ev.tr, data);   % add death mark
                    ev.save = 1;
                    ev.display = ['Track ' num2str(ev.tr) ' marked with death flag'];                    
                    
                case need_tr_s   % requires both a track and a selection as well
                    
                    [x,y] = gcaxy(gca);
                    ind = identifyobj(handles, x, y);
                    if (~ind)
                        ev.display = 'No object in frame';
                        continue
                    end
                    
                    switch inkey % switch on objects requiring track AND selection

                        case '5' % add object (or entire track) to existing track                            
                            tr = sch.tr(ev.tc,ind);
                            if (~tr)   % object not already assigned to a track
                                sch = sch.addtotrack(ev.tr, ev.tc, ind);
                                ev.display = ['Added new object to track ' num2str(ev.tr)];
                                ev.save = 1;
                            elseif (tr == ev.tr)
                                ev.display = ['Object already in track ' num2str(ev.tr)];
                            else   % join two tracks
                                if (~ev.direction)    % no d
                                    ev.display = 'No direction specified';
                                else
                                    sch = sch.jointrack(ev.tr, ev.tc-ev.direction, tr, ev.tc);
                                    ev.display = ['Appended track ' num2str(tr) ' to ' num2str(ev.tr)];
                                    ev.save = 1;
                                end
                            end                                   
                        case '\'   % split an object along a central line                                                                                   
                            obj = sch.obj(ev.tc,ind);            % tile number for the object to be divided
                            set(handles.figure1,'Name','Set dividing line segment');
                            figure(1);
                            h = imline(gca);   % handle to the line object
                            cs = h.getPosition;   % get the centroid positions
                            delete(h);        % now delete the line object
                            
                            xt1 = round(cs(1,1));      % x,y positions within the tile
                            yt1 = round(cs(1,2));
                            xt2 = round(cs(2,1));
                            yt2 = round(cs(2,2));
                            
                            [m1,n1,x1,y1] = tile2im_xy(handles.tile,xt1,yt1);
                            [m2,n2,x2,y2] = tile2im_xy(handles.tile,xt2,yt2);
                            
                            if ((m1 ~= m2)|(n1 ~= n2))   % end of line segments don't lie within the same tile
                                if ((m1 ~= obj.m)|(n1 ~= obj.n))     % line segment in different tile compared to object selected
                                    ev.display = 'Dividing line segment must lie completely within one tile';
                                    continue
                                end
                            end
                            objout = splitobj(obj,X,Y,x1,y1,x2,y2);
                            sch = sch.delobj(ev.tc,ind);   % remove old object
                            sch = sch.addobj(ev.tc,objout);   % add new objects                                                        
                            sch_old = handles.sch;        % manually save sch for undo
                            handles.sch = sch;              % need to update to allow new objects to be identified                            
                            ind2 = identifyobj(handles, x, y);
                            sch = sch.addtotrack(ev.tr, ev.tc, ind2);                               
                            handles.sch = sch;
                            ev.save = 0;                  % performed this save manually to enable undo 
                            ev.display = 'Object split';                            
                        case '/'   % split an object into two using the watershed algorithm                                                        
                            
                            objin = sch.obj(ev.tc,ind);   % this is the object to be split
                            objout = ws(objin,X,Y);   % watershed on the object                                 
                            sch = sch.delobj(ev.tc,ind);   % remove old object                            
                            sch = sch.addobj(ev.tc,objout);   % add new objects
                            sch_old = handles.sch;        % manually save sch for undo
                            handles.sch = sch;              % need to update to allow new objects to be identified                            
                            ind2 = identifyobj(handles, x, y);     % should now be unable to load deleted object
                            sch = sch.addtotrack(ev.tr, ev.tc, ind2);
                            handles.sch = sch;
                            ev.save = 0;                  % performed this save manually to enable undo                            
                            ev.display = 'Performed watershed';
                            
                        case 'z'  % specify child of given track                                                                                                                   
                            if (~sch.tr(ev.tc,ind))   % child object found, but not yet assigned to track
                                sch = sch.addtotrack(0,ev.tc,ind);
                                ev.save = 1;
                            else                                                                                    
                                tr_c = sch.tr(ev.tc,ind);                                
                                if (tr_c == ev.tr)  % child object cannot be the same as the parent object                                    % object            
                                    ev.display = 'Child cannot be made the same as parent track';
                                    continue
                                end                                   
                                [sch, added] = sch.addchild(ev.tr, tr_c);                                
                                ev.save = 1;
                                ev.display = ['Track ' num2str(ev.tr) ' children: ' num2str(sch.children(ev.tr))];                                
%                                 if added
%                                     ev.display = ['Track ' num2str(tr_c) ' added as child of track ' num2str(ev.tr)];
%                                     ev.save = 1;
%                                 else
%                                     ev.display = 'Track already has two children';
%                                 end
                            end
                    end
            end
    end
    
end
set(handles.figure1,'Name','Stack Viewer')
set(gcf,'WindowStyle','normal');
guidata(hObject,handles);
return
end

function [m,n,x,y] = tile2im_xy(tile, xnow, ynow)

Ym = length(tile.closem);
Xm = length(tile.closen);

ynow = min(ynow, Ym);
xnow = min(xnow, Xm);

m = tile.closem(ynow);
n = tile.closen(xnow);
x = tile.inx(ynow,xnow,m,n);
y = tile.iny(ynow,xnow,m,n);

end


function indpoint  = identifyobj(handles, xnow, ynow)

%% get location of the pixel within the image tile
tc = round(get(handles.planecommon,'Value'));   % the index for this common frame that is currently displayed
sch = handles.sch;

[m,n,x,y] = tile2im_xy(handles.tile,xnow,ynow);

obj = sch.timepoints(tc).obj;

%m = handles.tile.closem(ynow);   % the row number for the image to be segmented
%n = handles.tile.closen(xnow);   % the column number for the image to be segmented
%x = handles.tile.inx(ynow,xnow,m,n);   % x,y indices for the pixel location within the tile
%y = handles.tile.iny(ynow,xnow,m,n);
ms = [obj.m];   % find all the objects in the current frame
ns = [obj.n];
trnos = [obj.trno];
ind = intersect(find(ms==m), find(ns==n));
ind = intersect(ind, find(trnos > -1));   % exclude the killed cells

%% restrict selection to within gate, if gate is selected
selectedgate = (get(handles.gates_select,'Value') - 1);
if (selectedgate)
    cellgates = [obj.gate];    
    ind = intersect(ind, find(cellgates==selectedgate));
end
   
if (~isempty(ind));    
    xs = [obj(ind).x];
    ys = [obj(ind).y];
    trnos = [obj(ind).trno];    
    ds = (x-xs).^2 + (y-ys).^2;   % distance squared from the cursor
    imin = find(ds == min(ds));   % index of the object that has minimum distance from the cursor
    imin = imin(1);    
    indpoint = ind(imin);   % this is the index that is pointed to
else
    indpoint = 0;    % no cells were present in the field of view
end  
end


function loadschnitz_Callback(hObject, eventdata, handles)
[filename dirname] = uigetfile('sch.mat','Choose Schnitz file...');
load([dirname '/' filename]);
sch = Schnitz2(timepoints, tracks);
handles.sch = sch;
guidata(hObject,handles);
end

function calculate_Callback(hObject, eventdata, handles)
% lets first try to find a way to fill in all the blank data fields
sch = handles.sch;
imdir = handles.acq.imdir;
tcommon = handles.acq.tcommon;

for t = 1:sch.T   % loop through all timepoints    
    imall = [];   
    obj = sch.obj(t);   % retrieve object structure for this timepoint        
    
    if (isempty(obj))
        continue
    end
    trs = [obj.trno];    % retrieve the list of all track numbers for the objects
    for i = 1:length(obj)        
        %if (trs(i) <= 0) % object not assigned to a track or deleted
        %    continue        
        %elseif (~sch.approved(trs(i))) % object not yet approved
        %    continue
        if (~isempty(obj(i).data)) % object already contains data
            continue
        end
            
        % calculate the data
        if isempty(imall)  % load timelapse image in not present already
            imall = loadone(imdir,tcommon(t));
        end
        data = getdata(imall, obj(i), handles);
        sch = sch.addobjdata(t, i, data);    
    end
end
handles.sch = sch;
guidata(hObject, handles);

% modify this to obtain data for all approved cells
% then track them into lineages
sch = handles.sch;
imdir = handles.acq.imdir;
tcommon = handles.acq.tcommon;
thr = handles.acq.tr(tcommon);

% -1 check for self-referntial parent child relationships
% can remove this later in the next version of the tracking
for i = 1:length(sch.Tr)
    children = sch.tracks(i).children;
    if sum(ismember(i,children));
        sch.tracks(i).children = [];
        fprintf(['Track ' num2str(i) ' found to be self-referential.  Removing children.\n']);
    end
end

% 0 check graph to see that all parent/child relationships hold 
sch = sch.checklins();
% 1 approve all trajectories that have a parent/child
sch = sch.approveall();   
% 2 get trajectory data for all approved trajectories


if (~handles.tracked)
    for t = 1:sch.T    % loop through all the timepoints
        %waitbar(t/sch.T, h, ' Processing timepoints...');
        t
        imall = [];          % start with empty image
        obj = sch.obj(t);    % get all objects


        
        for i = 1:length(obj)    
           if ~sch.approved(obj(i).trno)  % skip object if not approved
               continue                      
           elseif isempty(imall)  % load timelapse image
                imall = loadone(imdir,tcommon(t));   
           end       
           data = getdata(imall, obj(i), handles);
           sch = sch.addobjdata(t, i, data);
        end
    end
end
handles.sch = sch;
guidata(hObject,handles);  % save the information first

% process all approved timepoints
alltracks = find([sch.approved]);
L = length(alltracks);   % the initial number of tracks to be processed 
i = 1;         % generate all the tracks
while (~isempty(alltracks))  
    
    1 - length(alltracks)/L    
    %waitbar(1 - length(alltracks)/L, h, 'Tracking lineages...');
    trs = sch.getlintr(alltracks(1));   
    alltracks(1)
    lin(i).trs = trs; % get the tracks for each lineage    
    tracked = [trs.tr];    
    [tracked, inda, indb] = intersect(alltracks, tracked);
    
    if (isempty(intersect(inda, 1)))
        fprintf('cannot retrace graph\n');
        pause;
    end    
    i = i+1;
    i    
    alltracks
    alltracks(inda) = [];
end
%close(h)
handles.tracked = 1;

[filename dirname] = uiputfile('lin.mat','Save lineages file');
save([dirname '/' filename], 'lin', 'thr');
handles.sch = sch;
%3 somehow find a way to process the timepionts.
end

function data = getdata(imall, obj, handles);
    
    X = handles.acq.X;
    Y = handles.acq.Y;
    
    Xo = str2num(get(handles.xoffset,'String'));   % offset values
    Yo = str2num(get(handles.yoffset,'String'));
    
    % generate binary mask
    inds = sub2ind([Y X], obj.b(:,1), obj.b(:,2));   % object indices
    im1 = zeros(Y,X);   % binary image with the selected objected
    im1(inds) = 1;            
    imb = imfill(im1,'holes');   % fill all the holes in the image first
    
    % basic object information - do we need to call a program to get this
    % information?
    data.area = sum(imb(:));   % object area
    data.perim = sum(im1(:));   % object perimeter    
    data.x = obj.x + (X-Xo).*(obj.n-1);      % x position
    data.y = obj.y + (Y-Yo).*(obj.m-1);      % y position with tiling      
    % integrated intensity information
    data.erkbfp = sum(sum(imall(2).im(:,:,obj.m,obj.n).*uint16(imb)));
    data.h2bifp = sum(sum(imall(3).im(:,:,obj.m,obj.n).*uint16(imb)));
    data.nfatruby = sum(sum(imall(4).im(:,:,obj.m,obj.n).*uint16(imb)));
end

function objout = newobj(imb, m, n)    % takes new object specified by binary image imb, and generates new object 

[b,iml] = bwboundaries(imb);  % make the boundary image using four connectivity
s = regionprops(iml,'Centroid');
objout.m = m;
objout.n = n;
objout.x = s.Centroid(1);
objout.y = s.Centroid(2);
objout.gate = 0;
objout.b = b{1};
objout.trno = 0;   
objout.data = [];
end

function oh = showobj(obj, tile, color, width, style)
xmin = tile.xmins(obj.m, obj.n);
ymin = tile.ymins(obj.m, obj.n);
x = obj.x+xmin-1;
y = obj.y+ymin-1;
oh = plot(obj.b(:,2)+xmin-1, obj.b(:,1)+ymin-1, 'Color', color,'LineWidth',width,'LineStyle',style);
end

function eraseline(tile)
figure(1);
allobj = findobj(gca,'Type','Line');
notborder = setdiff(allobj, [tile.bvh tile.bih]);
delete(notborder);  % delete all the object boundaries in the image
end

function images = loadone(imdir,i)   % load the timelapse image at frame i
filename = [imdir 'imgf_' num2strn(i,4) '.mat'];
load(filename);
if ~exist('images')
    error('Invalid image files...');
end
end
% --- Executes on button press in showborder.
function showborder_Callback(hObject, eventdata, handles)

tile = handles.tile;
if get(hObject, 'Value')
    tile = makeborder(tile);
else   % delete the border
    bvh = tile.bvh;
    bih = tile.bih;
    delete(bvh);
    delete(bih);
    tile.bvh = 0;
    tile.bih = 0;
end
handles.tile = tile;   % save all the tiling information back to the handles structure
guidata(hObject, handles);
end
function showmissing_Callback(hObject, eventdata, handles)
handles = refreshfig(handles,1);
guidata(hObject, handles);
end

% --- Executes on button press in showtime.
function showtime_Callback(hObject, eventdata, handles)
refreshfig(handles,0);
guidata(hObject,handles);
end
% --- Executes on button press in makemovie.
function makemovie_Callback(hObject, eventdata, handles)

T = handles.acq.T;   % all the planes in the movie
tnow = str2num(get(handles.plane_current,'String'));  % store the current frame of the movie so that we can return to it later

frames = inputdlg('Enter frames for movie:');
frames = eval(frames{1});
if (min(frames)<1) | (max(frames)>T)
    error('Movie frames out of range...');
end
fps = inputdlg('Enter frames per second form movie');
fps = str2num(fps{1});

[moviename, pathname] = uiputfile('*.avi', 'Enter name of output movie AVI file:');
  
    %Prepare the new file.
    vidObj = VideoWriter([pathname '/' moviename], 'Uncompressed AVI');
    vidObj.FrameRate = fps;
    open(vidObj);
    figure(1);
    for i = 1:length(frames)
        f = frames(i);
        set(handles.plane_slider,'Value',f);
        set(handles.plane_current,'String',num2str(f));
        set(handles.planecommon,'Value',f);
        handles = refreshfig(handles,1);  % reload every time        
        writeVideo(vidObj,getframe(gcf));
    end
    close(vidObj);

% set the plane back to the previous plane before the program started
set(handles.plane_slider,'Value',tnow);
set(handles.plane_current,'String',num2str(tnow));
handles = refreshfig(handles,0);
guidata(hObject, handles);
end

function tr_set(handles, tr, rgb);
if (~tr)
    rgb = [0.5 0.5 0.5];   % set a gray tone if no object is selected
end
set(handles.tr,'String',num2str(tr));
set(handles.tr,'BackgroundColor',rgb);
end

function ontop_Callback(hObject, eventdata, handles)
fprintf('Hello World!\n');
setAlwaysOnTop(handles.figure1, get(hObject,'Value'));
end

function saveschnitz_Callback(hObject, eventdata, handles)
sch = handles.sch;
timepoints = sch.timepoints;
tracks = sch.tracks;
gatenames = get(handles.gates_select,'String');

[filename dirname] = uiputfile('sch.mat','Save tracks.mat file');
save([dirname '/' filename], 'timepoints','tracks','gatenames');
end

function objdata_Callback(hObject, eventdata, handles)
end

% --- Executes on selection change in gates_select.
function gates_select_Callback(hObject, eventdata, handles)
refreshfig(handles,0);
end

function gates_select_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function objout = ws(objin,X,Y)
objinds = sub2ind([Y X], objin.b(:,1), objin.b(:,2));               
im = zeros(Y,X);   % binary image with the selected objected
im(objinds) = 1;            

im = imfill(im,'holes');   % fill all the holes in the image first
im = -bwdist(~im);
L = watershed(im);
im2 = (L & im);  % this is the watershed segmented image

[b,iml] = bwboundaries(im2);  % make the boundary image
s = regionprops(iml,'Centroid');
for i = 1:length(s)
    objout(i).m = objin.m;
    objout(i).n = objin.n;
    objout(i).x = s(i).Centroid(1);
    objout(i).y = s(i).Centroid(2);
    objout(i).gate = objin.gate;
    objout(i).b = b{i};
    objout(i).trno = 0;  
    objout(i).data = [];
end
end

function objout = splitobj(objin,X,Y,x1,y1,x2,y2);
% takes in object OBJIN, and coordinates X1,Y1,X2,Y2 that define a line
% then splits OBJ into multiple objects around the dividing line

r = 0:0.002:1;

xp = round(x1 + (x2-x1).*r + (rand-0.5).*3);   % parametrized version of them line
yp = round(y1 + (y2-y1).*r + (rand-0.5).*3);   % parametrized version of the line


inds_line = sub2ind([Y X], yp, xp);
inds_obj = sub2ind([Y X], objin.b(:,1), objin.b(:,2));

im = zeros(Y,X);
im(inds_obj) = 1;   % make dividing line;
im = imfill(im,'holes');   % fill holds;
im(inds_line) = 0;

[b,iml] = bwboundaries(im);  % make the boundary image using four connectivity
s = regionprops(iml,'Centroid');

for i = 1:length(s)
    objout(i).m = objin.m;
    objout(i).n = objin.n;
    objout(i).x = s(i).Centroid(1);
    objout(i).y = s(i).Centroid(2);
    objout(i).gate = objin.gate;
    objout(i).b = b{i};
    objout(i).trno = 0;   
    objout(i).data = [];
end
end


% --- Executes on button press in purge.
function purge_Callback(hObject, eventdata, handles)

    T = handles.acq.T;    % the total number of frames in the stack            
    for i = 1:T
        handles.buffer(i).imall=[];
    end    
    set(handles.figure1,'Name','buffer purged');   % display the text to be displayed
    guidata(hObject, handles);  % save
end

% hObject    handle to purge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
