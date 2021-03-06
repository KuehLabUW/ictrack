function varargout = plot_multicolor(varargin)
% PLOT_MULTICOLOR M-file for plot_multicolor.fig
%      PLOT_MULTICOLOR, by itself, creates a new PLOT_MULTICOLOR or raises the existing
%      singleton*.
%
%      H = PLOT_MULTICOLOR returns the handle to a new PLOT_MULTICOLOR or the handle to
%      the existing singleton*.
%
%      PLOT_MULTICOLOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLOT_MULTICOLOR.M with the given input arguments.
%
%      PLOT_MULTICOLOR('Property','Value',...) creates a new PLOT_MULTICOLOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before plot_multicolor_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to plot_multicolor_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help plot_multicolor

% Last Modified by GUIDE v2.5 28-Oct-2011 13:14:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @plot_multicolor_OpeningFcn, ...
                   'gui_OutputFcn',  @plot_multicolor_OutputFcn, ...
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


% --- Executes just before plot_multicolor is made visible.
function plot_multicolor_OpeningFcn(hObject, eventdata, handles, varargin)

if (nargin ~= 6)
    error('USAGE:  plot_multicolor(DATA,X,Y)');
end

data = varargin{1};
X = varargin{2};
Y = varargin{3};

handles.data = data;
handles.X = X;
handles.Y = Y;

refresh(handles);

% Choose default command line output for plot_multicolor
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


function refresh(handles);

data = handles.data;
X = handles.X;
Y = handles.Y;

N = length(data)   % the number of data points
imall = ones(length(Y)-1, length(X)-1, 3);    % make the color image, initially white

for i = 1:N
    Z = data(i).Z;
    Z = Z./max(Z(:));   % normalize maximum to unity    
    r = str2num(get(handles.(['r' num2str(i)]),'String'));
    g = str2num(get(handles.(['g' num2str(i)]),'String'));
    b = str2num(get(handles.(['b' num2str(i)]),'String'));    
    imone(:,:,1) = (g+b).*Z;
    imone(:,:,2) = (r+b).*Z;
    imone(:,:,3) = (r+g).*Z;
    %figure(i); imshow(imone)
    imall = imall - imone;    
end

figure(12);
imshow(imall,'XData',X,'YData',Y);
axis xy; axis on; axis normal;

return





% --- Outputs from this function are returned to the command line.
function varargout = plot_multicolor_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function r1_Callback(hObject, eventdata, handles)
refresh(handles);


function r1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function g1_Callback(hObject, eventdata, handles)
refresh(handles);


function g1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function b1_Callback(hObject, eventdata, handles)
refresh(handles);

function b1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function r2_Callback(hObject, eventdata, handles)
refresh(handles);

function r2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function g2_Callback(hObject, eventdata, handles)
refresh(handles);

function g2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function b2_Callback(hObject, eventdata, handles)
refresh(handles);

function b2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function r3_Callback(hObject, eventdata, handles)
refresh(handles);

function r3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function g3_Callback(hObject, eventdata, handles)
refresh(handles);


function g3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function b3_Callback(hObject, eventdata, handles)
refresh(handles);

function b3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function refresh_Callback(hObject, eventdata, handles)
refresh(handles)
