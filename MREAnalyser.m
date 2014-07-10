function varargout = MREAnalyser(varargin)
% MREANALYSER MATLAB code for MREAnalyser.fig
%      MREANALYSER, by itself, creates a new MREANALYSER or raises the existing
%      singleton*.
%
%      H = MREANALYSER returns the handle to a new MREANALYSER or the handle to
%      the existing singleton*.
%
%      MREANALYSER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MREANALYSER.M with the given input arguments.
%
%      MREANALYSER('Property','Value',...) creates a new MREANALYSER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MREAnalyser_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MREAnalyser_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MREAnalyser

% Last Modified by GUIDE v2.5 10-Jul-2014 12:17:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MREAnalyser_OpeningFcn, ...
                   'gui_OutputFcn',  @MREAnalyser_OutputFcn, ...
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


% --- Executes just before MREAnalyser is made visible.
function MREAnalyser_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MREAnalyser (see VARARGIN)

% Choose default command line output for MREAnalyser
handles.output = hObject;

% Set internal variables
handles.phase = 1; % The active phase in the GUI
handles.slice = 1; % The active slice in the GUI
handles.image = cell(1);  % An array of images in case several are open
handles.curImg = 1;  % Index of the current image
handles.images{handles.curImg} = struct('phase',[],'mag',[]);
  % Each image has both magnitude and phase information
handles.activeImg = handles.images{handles.curImg}.phase;
  % The image that is being displayed on screen

% Update handles structure
guidata(hObject, handles);

% Prompt for file
status(handles,'Open an image file')

% UIWAIT makes MREAnalyser wait for user response (see UIRESUME)
% uiwait(handles.figure1);
colormap gray;


% --- Outputs from this function are returned to the command line.
function varargout = MREAnalyser_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function sliceSelect_Callback(hObject, eventdata, handles)
% hObject    handle to sliceSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.slice = round(get(hObject,'Value'));
set(hObject,'Value',handles.slice);
guidata(hObject,handles);
redraw(handles);

% --- Executes during object creation, after setting all properties.
function sliceSelect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliceSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in deAliasBtn.
function deAliasBtn_Callback(hObject, eventdata, handles)
% hObject    handle to deAliasBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.statusBar,'String','dealiasing');

% --- Executes on slider movement.
function phaseSelect_Callback(hObject, eventdata, handles)
% hObject    handle to phaseSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.phase = round(get(hObject,'Value'));
set(hObject,'Value',handles.phase);
guidata(hObject,handles);
redraw(handles);




% --- Executes during object creation, after setting all properties.
function phaseSelect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phaseSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --------------------------------------------------------------------
function fileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to fileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function openMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to openMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
openFile(handles)


% --- Executes during object creation, after setting all properties.
function sliceWindow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliceWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate sliceWindow
axis square;

% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliceWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate sliceWindow

function redraw(handles)
s = handles.slice;
p = handles.phase;
c = get(handles.cAxis,'Value');
switch get(handles.magToggle,'SelectedObject')
    case handles.magBtn
        clim = [0 c];
    case handles.phaseBtn
        clim = [-c c];
end
imagesc(handles.activeImg(:,:,s,p),...
        'Parent',handles.sliceWindow);
    axis square;
    caxis(clim);


% --- Executes on mouse press over axes background.
function sliceWindow_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to sliceWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in playMovie.
function playMovie_Callback(hObject, eventdata, handles)
% hObject    handle to playMovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
makeSliceMovie(handles,false)


% --- Executes on key press with focus on phaseSelect and none of its controls.
function phaseSelect_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to phaseSelect (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
%disp(eventdata)


% --- Executes on slider movement.
function cAxis_Callback(hObject, eventdata, handles)
% hObject    handle to cAxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
scale = get(hObject,'Value');
switch get(handles.magToggle,'SelectedObject')
    case handles.magBtn
        set(handles.sliceWindow,'CLim',[0,scale]);
    case handles.phaseBtn
        set(handles.sliceWindow,'CLim',[-scale,scale]);
end

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function cAxis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cAxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --------------------------------------------------------------------
function openFile_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to openFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
openFile(handles)

function openFile(handles)
f = handles.sliceWindow;
i = handles.curImg;
status(handles,'Loading...');
[mag, phase] = getMREimages(1:30,false,f);
status(handles,'Image loaded')
handles.images{i}.mag = mag;
handles.images{i}.phase = phase;
handles.activeImg = handles.images{i}.phase;
p = handles.phaseSelect;
s = handles.sliceSelect;
dims = size(handles.activeImg);
numSlices = dims(3); numPhases = dims(4);
set(s,'Min',1), set(s,'Max',numSlices), set(s, 'Value',1);
set(p,'Min',1), set(p,'Max',numPhases), set(p, 'Value',1);
set(s, 'SliderStep', [1 1]/(numSlices-1));
set(p, 'SliderStep', [1 1]/(numPhases-1));

redraw(handles);
guidata(hObject,handles);
set(handles.magToggle,'SelectedObject',handles.phaseBtn);
updateCAxis(handles);


% --------------------------------------------------------------------
function saveMovie_Callback(hObject, eventdata, handles)
% hObject    handle to MovieMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
makeSliceMovie(handles,true)

function makeSliceMovie(handles,save)
colormap gray, 
slice = handles.slice;
for i = 1:get(handles.phaseSelect,'Max'),
    set(handles.phaseSelect,'Value',i);
    handles.phase = i;
    
    %guidata(handles.phaseSelect,handles);
    %imagesc(handles.activeImg(:,:,slice,i),...
      %'Parent',handles.sliceWindow);
      redraw(handles)
    if save
        M(i) = getframe;
    else
        pause(.2);
    end
end
if save
    fname = inputdlg('Filename to save movie:');
    fname = sprintf('%s.avi',fname{1});
    movie2avi(M,fname);
end


% --- Executes when selected object is changed in magToggle.
function magToggle_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in magToggle 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
i = handles.curImg;
switch hObject
    case handles.magBtn
        handles.activeImg = handles.images{i}.mag;
    case handles.phaseBtn
        handles.activeImg = handles.images{i}.phase;
end
guidata(hObject,handles)
updateCAxis(handles)
redraw(handles)



function status(handles,msg)
% Sets status bar text
set(handles.statusBar,'String',msg);

function updateCAxis(handles)
% Resets sliders for color
i = handles.curImg;
switch get(handles.magToggle,'SelectedObject')
    case handles.magBtn
        m = abs(handles.images{i}.mag);
    case handles.phaseBtn
        m = abs(handles.images{i}.phase);
end

while numel(m) > 1
    m = max(m);
end

set(handles.cAxis,'Value',0);
set(handles.cAxis,'Max',1.5*m);
set(handles.cAxis,'Value',0.75*m);
