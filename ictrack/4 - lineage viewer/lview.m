function varargout = lview(varargin)
% LVIEW M-file for lview.fig
%      LVIEW, by itself, creates a new LVIEW or raises the existing
%      singleton*.
%
%      H = LVIEW returns the handle to a new LVIEW or the handle to
%      the existing singleton*.
%
%      LVIEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LVIEW.M with the given input arguments.
%
%      LVIEW('Property','Value',...) creates a new LVIEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before lview_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to lview_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help lview

% Last Modified by GUIDE v2.5 02-Jun-2011 12:49:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @lview_OpeningFcn, ...
                   'gui_OutputFcn',  @lview_OutputFcn, ...
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

function lview_OpeningFcn(hObject, eventdata, handles, varargin)

addpath ../bin;
handles.output = hObject;
% load the lin.mat file
[imfile imdir] = uigetfile('*.mat', 'Point to lin.mat file...');
load([imdir imfile]);   
if ~exist('lin')
    error('Cannot find lin structure in file...');
end
% initialize the slider parameters
L = length(lin);
set(handles.linslider,'Value',1);
set(handles.linslider,'Min',1);
set(handles.linslider,'Max',L);
set(handles.linslider,'SliderStep',[1/(L-1) 5/(L-1)]);

% initialize the data fields class, if it doesn't already exist
if ~exist('df')
    df = DFS(fieldnames(lin(1).trs(1).data));
end
set(handles.fnameselect, 'String', ['none' df.fnames]);


% now we need to convert all the data for structure format to vector
% format
for i = 1:length(lin)
    trs = lin(i).trs;
    for j = 1:length(trs)
        data = trs(j).data;   % handles structure                
        if (isempty(data))
            fprintf(['no data found for lin ' num2str(i) ' track ' num2str(j) '.\n']);
            continue
        end
        for k = 1:length(df.fnames)   % loop through all the field names            
            if isfield(data, df.fnames{k})
                datanew.(df.fnames{k}) = [data.(df.fnames{k})];
            end
        end
        lin(i).trs(j).data = datanew;
    end
end
% the three variables to be saved and processed
handles.lin = lin;    % this is the raw lineage information 
handles.df = df;      % the list of fields derived for the data
handles.thr = thr;    % this is the time scale 

guidata(hObject, handles);

function varargout = lview_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function fnameselect_Callback(hObject, eventdata, handles)
    df = handles.df;
    fname = fnamefind(handles);    
    
    if strcmp(fname,'none');
        set(handles.displaymin,'Enable','inactive');
        set(handles.displaymax,'Enable','inactive');        
    else           
        set(handles.displaymin,'String',num2str(df.getmin(fname)));    
        set(handles.displaymax,'String',num2str(df.getmax(fname)));
        set(handles.displaymin,'Enable','on');
        set(handles.displaymax,'Enable','on');        
    end
    handles = refreshfig(handles);

function fnameselect_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function linslider_Callback(hObject, eventdata, handles)
i = round(get(hObject,'Value'));
set(hObject,'Value',i);   % ensure that the value is a round number
set(handles.linnum,'String',num2str(i));
handles = refreshfig(handles);

function linslider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function displaymin_Callback(hObject, eventdata, handles)
    val = str2num(get(hObject,'String'));
    df = handles.df;
    fname = fnamefind(handles);    
    handles = refreshfig(handles);
    handles.df = handles.df.setmin(fname,val);
    guidata(hObject,handles);

function displaymin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function displaymax_Callback(hObject, eventdata, handles)
    val = str2num(get(hObject,'String'));
    df = handles.df;
    fname = fnamefind(handles);
    handles.df = handles.df.setmax(fname,val);
    handles = refreshfig(handles);
    guidata(hObject,handles);

function displaymax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function makehistt_Callback(hObject, eventdata, handles)
    fname = fnamefind(handles);    
    if strcmp(fname,'none');
        set(handles.figure1,'String','No attribute selected for histogram');
        return
    end
    
    thr =  handles.thr;    % this is the elapsed time in hours
    lin = handles.lin;
    T = length(thr);       % this is the number of timepoints
    
    timepoints(T).data = [];   % initialize the empty structure
        
    for i = 1:length(lin);    % loop through all the lineages        
        i
        trs = lin(i).trs;     % loop through all tracks within a given lineage                
        for j = 1:length(trs)
            
            ts = trs(j).ts;  
            if isempty(ts)  % check if timepoints are empty
                continue
            end            
            datas = [trs(j).data.(fname)];            
            for k = 1:length(ts)
                timepoints(ts(k)).data =  [timepoints(ts(k)).data log10(datas(k))];
            end
        end
    end
    
    % calculate histogram over all possible timepoints
    alldata = [timepoints.data];    % this is the vector of all values over all timepoints
%     c = prctile(alldata, [5 50 95]);   % calculate median and the 10 and 90th percentiles    
%     clow = c(2) - 1.5*(c(2)-c(1));
%     chigh = c(2) + 1.5*(c(3)-c(2));    
    clow = 3.2; %log10(1000);
    chigh = 5.2; %log10(6*10^5);

    bins = clow:(chigh-clow)./100:chigh;   % bins for the histogram    
    N = hist(alldata,bins);
    figure(11); area(bins,N); 
    xlabel('PU.1-GFP levels');
    ylabel('frequency');
    title('PU.1-GFP level distribution over all times');
    
    % now bin measurements over time to help with sampling frequency
    tbinsize = 5;    
    tbins = floor(T/tbinsize);
    tbins = repmat(1:tbins, tbinsize, 1);
    tbins = tbins(:);    
    tbinned(tbins(end)).data = [];
        for i = 1:length(tbins)        
        tbinned(tbins(i)).data = [tbinned(tbins(i)).data timepoints(i).data];
    end           
    % now construct histogram for each binned timepoint    
    thist = zeros(length(bins), length(tbinned));    
    for i = 1:length(tbinned)
        thist(:,i) = (hist(tbinned(i).data, bins))';
    end
    
    thrb = thr((1:length(tbinned)).*tbinsize);    
    figure(12); surf(bins, thrb, thist');
    shading interp;
    axis off
    axis ij;
    view(2);
    keyboard
    %print('-djpeg', '-r90','-painters', jpgfilename)
    % save the 
    
    
    
    
function handles = refreshfig(handles)
    
    i = str2num(get(handles.linnum,'String'));
    lin = handles.lin;  
    df = handles.df;
    thr = handles.thr;        
    fname = fnamefind(handles);    
    if strcmp(fname,'none');        
        figure(1); cla;
        plotlin(thr, lin(i).trs);
    else
        figure(1); cla;
        plotlin(thr, lin(i).trs, fname, df.getcolor(fname), df.getmax(fname));                
        axis([0 100 -1.2 1.2]);         %hack
        figure(2); cla;
        plotlinflat(thr, lin(i).trs, fname, df.getcolor(fname), df.getmax(fname));        
        axis([0 100 0 0.3]);   %hack       
    end
    

function fname = fnamefind(handles)
        select = get(handles.fnameselect,'Value');
        fnames = get(handles.fnameselect,'String');
        fname = fnames{select};

function newfield_Callback(hObject, eventdata, handles)
fn = inputdlg('Enter output field name:');
fn = fn{1};
fnames = get(handles.fnameselect,'String');

for i = 1:length(fnames);
    fprintf(['{' num2str(i) '} - ' fnames{i} '\n']);
end
fprintf(['{t} - time\n']);

expr = inputdlg('Enter expression, enclose fields in braces:');
expr = expr{1};

% now need to replace expressions here with the dynamic field names

thr = handles.thr;

for i = 2:length(fnames)
    instr = ['{' num2str(i) '}'];
    outstr = ['trs(j).data.(fnames{' num2str(i) '})'];
    expr = regexprep(expr,instr,outstr);
end
% replace the time variable
expr = regexprep(expr,'{t}', 'thr(trs(j).ts)');
lin = handles.lin;
for i = 1:length(lin)        
    trs = lin(i).trs;
    for j = 1:length(trs)              
        if isempty(trs(j).ts)
            continue   % no data!
        elseif (length(trs(j).ts)==1)
            % don't process single data points
            continue            
        end        
        data = trs(j).data;                
        trs(j).data.(fn) = eval(expr);        
    end
    lin(i).trs = trs;
end
% assign new field name
fnames{end+1} = fn;
set(handles.fnameselect,'String',fnames);
handles.df = handles.df.addfield(fn);
handles.lin = lin;
guidata(hObject,handles);

function save_Callback(hObject, eventdata, handles)
lin = handles.lin;
thr = handles.thr;
df = handles.df;
[file path] = uiputfile('*.mat', 'Save lineage file...');
save([path file], 'lin', 'thr','df');

function button1_Callback(hObject, eventdata, handles)

lin = handles.lin;
df = handles.df;
thr = handles.thr;        
fname = fnamefind(handles);    

traces = [1 2 3 5 9 12 14 17 18 19 21 22 24 25 32 33 43 44 54 56 57 60 61 72];
count = 1;

for i = traces    
    if strcmp(fname,'none');
        fnumgen(1,3,6,count)
        plotlin(thr, lin(i).trs);
    else
        fnumgen(1,2,6, count)
        plotlin(thr, lin(i).trs, fname, df.getcolor(fname), df.getmax(fname));
        set(gca,'YTick',[]); title([num2str(i)],'FontSize',12);
        axis([0 100 -1.2 1.2]);         %hack
        fnumgen(101,2,6, count)
        plotlinflat(thr, lin(i).trs, fname, df.getcolor(fname), df.getmax(fname));
        set(gca,'YTick',[]); title([num2str(i)],'FontSize',12);
        axis([0 100 0 0.6]);   %hack
    end
    count = count+1;
end
 
function fnumgen(start,m,n,i)
    MN = m*n;   % total number of subplots per plot
    fignum = start + floor(i/MN);    
    snum = mod(i,MN);
    if (snum == 0)
        snum = MN;
        fignum = fignum-1;
    end
    figure(fignum);
    subplot(m,n,snum);    
    

    
    
    
    

