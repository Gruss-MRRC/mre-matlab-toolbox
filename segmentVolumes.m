function varargout = segmentVolumes(varargin)
% SEGMENTVOLUMES MATLAB code for segmentVolumes.fig
%      SEGMENTVOLUMES, by itself, creates a new SEGMENTVOLUMES or raises the existing
%      singleton*.
%
%      H = SEGMENTVOLUMES returns the handle to a new SEGMENTVOLUMES or the handle to
%      the existing singleton*.
%
%      SEGMENTVOLUMES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEGMENTVOLUMES.M with the given input arguments.
%
%      SEGMENTVOLUMES('Property','Value',...) creates a new SEGMENTVOLUMES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before segmentVolumes_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to segmentVolumes_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help segmentVolumes

% Last Modified by GUIDE v2.5 24-Jul-2014 17:18:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @segmentVolumes_OpeningFcn, ...
                   'gui_OutputFcn',  @segmentVolumes_OutputFcn, ...
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


% --- Executes just before segmentVolumes is made visible.
function segmentVolumes_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to segmentVolumes (see VARARGIN)

% Choose default command line output for segmentVolumes
handles.output = hObject;
handles.M = [];
pt = [1,1,1,1];
handles.pt = pt;
handles.volume = 1;
handles.cut = {[],[],[],[]};
handles.ptSelect =1;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes segmentVolumes wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = segmentVolumes_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in openButton.
function openButton_Callback(hObject, eventdata, handles)
% hObject    handle to openButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[M,~] = getRawMRE();
handles.M = M;

% Set volume selection slider properties based on number of volumes
% detected
dims = size(M);
set(handles.volumeSelect,'Max',dims(4));
set(handles.volumeSelect,'Min',1)
set(handles.volumeSelect,'SliderStep',[1 1]/dims(4));
pt = round(dims/2);
set(handles.volumeSelect,'Value',pt(4));
handles.pt = pt;
guidata(hObject, handles);
update_Axes(hObject,handles);


% --- Executes on button press in calculateButton.
function calculateButton_Callback(hObject, eventdata, handles)

M = handles.M;
d = size(M);
L = labelmatrix(bwconncomp(M(:,:,:,2)));
L(L(:)==0) = NaN; % set all background (0) values to NaN to avoid counting
% disp('<Volume Histogram>')
% group_histogram = histc(L(:),0:max(L(:)))
% [~,index] = max(group_histogram);
% disp('</Volume Histogram>')
Lslice = L(:,:,round(d(3)*2/3));
h = figure; 
colormap jet
imagesc(Lslice,[0 10]) % TODO don't hardcode the scale
coords = round(ginput(h));
val = Lslice(coords(2),coords(1));
volume = sum(L(:)==val);
close(h)
msgbox(sprintf('Volume of selected region: %g mL',volume/1000))
% bwselect();





% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% XZ plane
newpt = get(hObject,'CurrentPoint');
p = handles.figure1;
oldpt = handles.pt;
switch get(p,'SelectionType')
  case 'normal'%left mouse button click
    handles.pt = round([newpt(1),oldpt(2),newpt(3),oldpt(4)]);   
  case 'alt'%right mouse button click
    fprintf(2,'right click\n')
    i = handles.ptSelect;
    handles.cut{i} = round([newpt(1),oldpt(2),newpt(3),oldpt(4)]);
    if handles.ptSelect<4
      handles.ptSelect=handles.ptSelect+1; 
    else
      handles.ptSelect=1;
    end
end
guidata(hObject, handles);
update_Axes(hObject,handles);


% --- Executes on mouse press over axes background.
function axes2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% XY Plane
newpt = get(hObject,'CurrentPoint');
p = handles.figure1;
oldpt = handles.pt;
switch get(p,'SelectionType')
  case 'normal'%left mouse button click
    
    handles.pt = round([newpt(3),newpt(1),oldpt(3),oldpt(4)]);   
  case 'alt'%right mouse button click
    fprintf(2,'right click\n')
    i = handles.ptSelect;
    handles.cut{i} = round([newpt(3),newpt(1),oldpt(3),oldpt(4)]);
    if handles.ptSelect<4
      handles.ptSelect=handles.ptSelect+1; 
    else
      handles.ptSelect=1;
    end
end
guidata(hObject, handles);
update_Axes(hObject,handles);


% --- Executes on mouse press over axes background.
function axes3_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% YZ Plane, transposed so Z is up
newpt = get(hObject,'CurrentPoint');
p = handles.figure1;
oldpt = handles.pt;
switch get(p,'SelectionType')
  case 'normal'%left mouse button click
        handles.pt = round([oldpt(1),newpt(1),newpt(3),oldpt(4)]);   
  case 'alt'%right mouse button click
    fprintf(2,'right click\n')
    i = handles.ptSelect,4;
    handles.cut{i} = round([newpt(1),oldpt(2),newpt(3),oldpt(4)]);
    if handles.ptSelect<4
      handles.ptSelect=handles.ptSelect+1; 
    else
      handles.ptSelect=1;
    end
end
guidata(hObject, handles);
update_Axes(hObject,handles);


function update_Axes(hObject,handles)
f = [handles.axes1,handles.axes2,handles.axes3];
M = handles.M;
pt = handles.pt;
colormap gray
axis tight
% Check for active cut and set values to 0 along the cut line
c = handles.cut;
if ~isempty(c{1}) && ~isempty(c{2}) && ~isempty(c{3}) && ~isempty(c{4}) 
  M(min([c{1}(1),c{2}(1),c{3}(1),c{4}(1)]) : max([c{1}(1),c{2}(1),c{3}(1),c{4}(1)]),...
    min([c{1}(2),c{2}(2),c{3}(2),c{4}(2)]) : max([c{1}(2),c{2}(2),c{3}(2),c{4}(2)]),...
    min([c{1}(3),c{2}(3),c{3}(3),c{4}(3)]) : max([c{1}(3),c{2}(3),c{3}(3),c{4}(3)]),...
    c{1}(4)) = 0;
  handles.cut = {[],[],[],[]};
  handles.M = M;
  guidata(hObject,handles)
end
imagesc(transpose(squeeze(M(:,pt(2),:,pt(4)))),'Parent',f(1),'HitTest','off');
imagesc(squeeze(M(:,:,pt(3),pt(4))),'Parent',f(2),'HitTest','off');
imagesc(transpose(squeeze(M(pt(1),:,:,pt(4)))),'Parent',f(3),'HitTest','off');
for i=1:3, axis(f(i),'tight'); end
guidata(hObject,handles)


% --- Executes on slider movement.
function volumeSelect_Callback(hObject, eventdata, handles)
% hObject    handle to volumeSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
v = round(get(hObject,'Value'));
handles.pt(4) = v;
guidata(hObject,handles)
update_Axes(hObject,handles)


% --- Executes during object creation, after setting all properties.
function volumeSelect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to volumeSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
