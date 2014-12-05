function varargout = MREBulkMotionCalculator(varargin)
% Calculate bulk motion amplitude of a cerebral MRE test sequence
% GUI to unwrap MRE test sequence and determine of bulk motion amplitude is
% sufficient to carry out the full sequence. 'MRE Viewer' is a version of
% this program with more functionality.
%
% Inputs:
% 4D DICOM file with .dicom extension (loaded using the 'Open File' button in the
% toolbar)
%
% Outputs:
% The scaled motion amplitude is printed to the GUI window
%
% Use the 'Save' button to store the processed image as multidimensional
% matrices in Matlab's .mat format. You can load this file later in Matlab
% from the file picker dialog, or pull the variables into the Matlab
% workspace using the 'load' command (useful for creating your own movies
% or image sequences).
%
% See also: GUIDE, GUIDATA, GUIHANDLES, FASTUNWRAP, MRI2MAT, MRESLICE2AVI
% 
% Authors
% -------
% Alex Krolick <amk283@cornell.edu>
% Mark Wagshul <mark.wagshul@einstein.yu.edu>

% Last Modified by GUIDE v2.5 25-Sep-2014 12:35:17
%
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MREBulkMotionCalculator_OpeningFcn, ...
                   'gui_OutputFcn',  @MREBulkMotionCalculator_OutputFcn, ...
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

% --- Executes just before MREBulkMotionCalculator is made visible.
function MREBulkMotionCalculator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MREBulkMotionCalculator (see VARARGIN)

% Choose default command line output for MREBulkMotionCalculator
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
statusMsg(handles,'Open an image file')

% UIWAIT makes MREBulkMotionCalculator wait for user response (see UIRESUME)
% uiwait(handles.figure1);
colormap gray;


% --- Outputs from this function are returned to the command line.
function varargout = MREBulkMotionCalculator_OutputFcn(hObject, eventdata, handles) 
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
openFile(hObject,handles)


% --- Executes during object creation, after setting all properties.
function sliceWindow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliceWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate sliceWindow
axis equal;

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
    axis equal off;
    box off;
    
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


function openFile_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to openFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
openFile(hObject,handles)


function openFile(hObject,handles)
% Open a MRE image file and pass it to the appropriate handler function

statusMsg(handles,'Opening...');
[f p] = uigetfile('*.dicom')
[mag,phase,info] = loadimage(p,f);
handles.header = info;
i = handles.curImg;
nSlices = size(mag,3);
nPhases = size(mag,4);
statusMsg(handles,'Opening image...');
handles.images{i}.mag = mag;
handles.images{i}.phase = phase;
handles.activeImg = phase;
set(handles.magToggle,'Visible','on');
statusMsg(handles,'Setting up workspace...');
p = handles.phaseSelect;
s = handles.sliceSelect;
handles.phase = 1;
handles.slice = 1;
set(s,'Min',1), set(s,'Max',nSlices), set(s, 'Value',handles.slice);
set(p,'Min',1), set(p,'Max',nPhases), set(p, 'Value',handles.phase);
set(s, 'SliderStep', [1 1]/max((nSlices-1),1));
set(p, 'SliderStep', [1 1]/max((nPhases-1),1));
guidata(hObject,handles);
set(handles.magToggle,'SelectedObject',handles.phaseBtn);
updateCAxis(handles);
redraw(handles);
statusMsg(handles,'Image loaded')


 function [mag,phase,info] = loadimage(p,f)
% Take a 4D DICOM file and output 5D matrices for phase and magnitude by
% shuffling the 4th dimension. Assume that the input format is [x,y,t],
% where the time axis is filled according to the sequence
% [(mag,slice1,dir1,t1),(mag,slice1,dir1,t2)...(phase,sliceN,dirN,tN)].
% Decompose into a matrix for phase and a matrix for magnitude, which
% should each be of the form [x,y,z,phase,direction], where phase is a
% temporal dimension. The difference between READNIFTI5D and READDICOM5D is
% that the DICOM header contains information about the number of
% dimensions, but the NIFTI header does not. The number of directions is
% assumed to be 4 in READNIFTI5D. Image orientations may also be different
% if using DICOM.
%
% See also getMRESinkus.m, readNIFTI5D, readDICOM5D_Mask
  
  % Open file
  im = dicomread([p f]);
  im = double(squeeze(im)); % squeeze gets rid of singleton dimensions
  
  % Metadata
  if exist('/gmrrc/mrbin/dcmdump','file') && isunix()
    [~,result] = system(['dcmdump ' p f...
      ' | grep NumberOfTemporalPositions'...
      '| awk ''{print $3}'' | head -n 1 | sed ''s/\[//g'' | sed ''s/]//g''']);
    nDirs = str2num(result);
    [~,result] = system(['dcmdump ' p f ...
      ' | grep "(2001,1018)" | awk ''{print $3}''']);
    nSlices = str2num(result);
    info = struct();
  else
    info = dicominfo([p f]);
    nDirs   = info.Private_2001_1081;
    nSlices = info.Private_2001_1018;
  end
  nX = size(im,1); 
  nY = size(im,2);
  nVolumes = size(im,3);
  nPhases = nVolumes / 2 / nDirs / nSlices;
  
  % Fill 5D matrices with volumes from 3D matrix
  k1=1;
  mag = zeros(nX,nY,nSlices,nPhases,nDirs);
  phase = mag;
  for slice = 1:nSlices,
    for dir = 1:nDirs,
      for ph = 1:nPhases,
        mag(:,:,slice,ph,dir) = im(:,:,k1);
        phase(:,:,slice,ph,dir) = im(:,:,nVolumes/2+k1);
        k1 = k1 + 1;
       end
    end
  end
  phase = double(phase)*pi/2048-pi;


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
        if length(size(handles.images{i}.phase))==5
            set(handles.dirBtnGrp,'Visible','off');
        end
    case handles.phaseBtn
        handles.activeImg = handles.images{i}.phase;
        if length(size(handles.images{i}.phase))==5
            set(handles.dirBtnGrp,'Visible','on');
        end
end
guidata(hObject,handles)
updateCAxis(handles)
redraw(handles)

% --------------------------------------------------------------------
function statusMsg(handles,msg)
% Sets status bar text
set(handles.statusBar,'String',msg);
pause(.01)

% --------------------------------------------------------------------
function updateCAxis(handles)
% Resets sliders for color
i = handles.curImg;
switch get(handles.magToggle,'SelectedObject')
    case handles.magBtn
        m = abs(handles.images{i}.mag);
    case handles.phaseBtn
        m = abs(handles.images{i}.phase);
end

m=max(m(:));

set(handles.cAxis,'Value',0);
set(handles.cAxis,'Max',1.5*m);
set(handles.cAxis,'Value',0.5*m);


% --- Executes on button press in helpBtn.
function helpBtn_Callback(hObject, eventdata, handles)
% hObject    handle to helpBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox(['Contact Mark Wagshul <mark.wagshul@einstein.yu.edu> ',...
  'or Alex Krolick <amk283@cornell.edu>'],'Help')


% --- Executes when selected object is changed in dirBtnGrp.
function dirBtnGrp_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in dirBtnGrp 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
i = handles.curImg;
switch hObject
    case handles.xBtn
        handles.activeImg = handles.images{i}.phase(:,:,:,:,1);
    case handles.yBtn
        handles.activeImg = handles.images{i}.phase(:,:,:,:,2);
    case handles.zBtn
        handles.activeImg = handles.images{i}.phase(:,:,:,:,3);
    case handles.bBtn
        handles.activeImg = handles.images{i}.phase(:,:,:,:,4);
end
guidata(hObject,handles)
updateCAxis(handles)
redraw(handles)


% --- Executes on button press in unwrapBtn.
function unwrapBtn_Callback(hObject, eventdata, handles)
% hObject    handle to unwrapBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
statusMsg(handles,'Phase unwrapping...')
i = handles.curImg;
im_mag = handles.images{i}.mag;
im_phase = handles.images{i}.phase;

tic
im_unwrap = fastUnwrap(im_mag,im_phase);
runtime=toc;
    
statusMsg(handles,sprintf('Unwrap took %.2f s',runtime))
handles.images{i}.phase = squeeze(im_unwrap);
handles.images{i}.mag = squeeze(im_mag);
handles.activeImg = handles.images{i}.phase(:,:,:,:,1);
updateCAxis(handles);
guidata(hObject,handles);
redraw(handles);


function saveBtn_ClickedCallback(hObject, eventdata, handles)
% Save to disk
header = handles.header;
matfile   = [header.path strtok(header.filename,'.') '.mat'];
M = handles.images{handles.curImg}.mag;
P = handles.images{handles.curImg}.phase;
save(matfile, 'header', 'M', 'P');
statusMsg(handles,sprintf('Saved %s',matfile))
guidata(hObject,handles);

function bulkmotionBtn_Callback(hObject, eventdata, handles)
% Call bulk motion calculation function & report result

% Conversion factors
% microns per radian of phase
scale = str2num(get(handles.conversionScaleBox,'String'));

% Set up
P = handles.images{handles.curImg}.phase;

% Call calculation function
amplitude = bulkmotion(P,scale,1,false);
cutoffamplitude = 40; % microns

% Report back to user
message = sprintf(['Bulk motion amplitude: %.2g microns\n'...
                   'Minumum needed for MRE: %.2g microns'],...
                   amplitude,cutoffamplitude);
statusMsg(handles,message)
guidata(hObject,handles);


function conversionScaleBox_Callback(hObject, eventdata, handles)


function conversionScaleBox_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
