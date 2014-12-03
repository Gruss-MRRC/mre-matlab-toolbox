function varargout = MREView(varargin)
% Preview and phase unwrap MRE sequences.
% 
% Use the 'Save' button to store the processed image as multidimensional
% matrices in Matlab's .mat format. You can load this file later in MREView
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

% Last Modified by GUIDE v2.5 24-Sep-2014 18:15:55
%
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MREView_OpeningFcn, ...
                   'gui_OutputFcn',  @MREView_OutputFcn, ...
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

% --- Executes just before MREView is made visible.
function MREView_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MREView (see VARARGIN)

% Choose default command line output for MREView
handles.output = hObject;

% Set internal variables
handles.phase = 1; % The active phase in the GUI
handles.slice = 1; % The active slice in the GUI
handles.direction = 1; % The active plane in the 5th dimension
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

% UIWAIT makes MREView wait for user response (see UIRESUME)
% uiwait(handles.figure1);
colormap gray;


% --- Outputs from this function are returned to the command line.
function varargout = MREView_OutputFcn(hObject, eventdata, handles) 
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

[mag,phase,info] = mri2mat();
handles.header = info;
i = handles.curImg;
nSlices = size(mag,3);
nPhases = size(mag,4);
statusMsg(handles,'Opening image...');
switch length(size(mag))
  case 4 % MRE image
    handles.images{i}.mag = mag;
    handles.images{i}.phase = phase;
    handles.activeImg = phase;
    set(handles.dirBtnGrp,'Visible','off');  
    set(handles.magToggle,'Visible','on');
  case 5 % Motion-sensitized MRE image
    handles.images{i}.mag = mag;
    handles.images{i}.phase = phase;
    handles.activeImg = phase;
    set(handles.dirBtnGrp,'Visible','on');
    set(handles.magToggle,'Visible','on');
end
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
      d = 1;
    case handles.yBtn
      d = 2;  
    case handles.zBtn
      d = 3;
    case handles.bBtn
      d = 4;  
end
handles.direction = d;
handles.activeImg = handles.images{i}.phase(:,:,:,:,d);
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
dims = size(handles.images{i}.phase);
nSlices= dims(3); nPhases= dims(4);
im_mag = handles.images{i}.mag;
im_phase = handles.images{i}.phase;
im_unwrap = im_phase;

method = questdlg('Use Fast (nearest-neighbor), Slow (quality-guided), or Constantini process to unwrap:', ...
 'Unwrap Algorithm', ...
 'Fast','Slow','Constantini','Fast');

%-------------------------5D-------------------------%
im1a_ph = handles.images{i}.phase; % initialize
if length(dims) == 5,
  switch method
    case 'Slow'
      w=waitbar(0,'Progress');
      tic
      im1a_ph = im_phase;
      for dir = 1:dims(5);
         for slice = 1:nSlices,
              for phase = 1:nPhases,
                  im1a_ph(:,:,slice,phase,dir) = QualityGuidedUnwrap2D_r1(...
                                    squeeze(im_mag(:,:,slice,phase,dir)),...
                                    squeeze(im_phase(:,:,slice,phase,dir)));
                 stat = sprintf('Slice: %g/%g, Phase: %g/%g, Dir: %g',...
                     slice,nSlices,phase,nPhases,dir); 
                 waitbar((dir-1+(slice-1+phase/nPhases)/nSlices)/dims(5),w,stat)
              end
          end
      end
      delete(w)
      runtime=toc;
    case 'Fast'
      tic
      im1a_ph = fastUnwrap(im_mag,im_phase);
      runtime=toc;   
    case 'Constantini'
      tic
      for dir = 1:dims(5);
        im1a_ph(:,:,:,:,dir) = cunwrap4D(im_phase(:,:,:,:,dir));
      end
      runtime=toc;
  end
     % Subtract phase due to background field 
    im_ph_P = im1a_ph(:,:,:,:,1)-im1a_ph(:,:,:,:,4); % Phase direction
    im_ph_M = im1a_ph(:,:,:,:,2)-im1a_ph(:,:,:,:,4); % Magnitude direction
    im_ph_S = im1a_ph(:,:,:,:,3)-im1a_ph(:,:,:,:,4); % Slice direction

    % Time (phase) average
    mn_P = mean(im_ph_P,4);
    mn_M = mean(im_ph_M,4);
    mn_S = mean(im_ph_S,4);

    % Direction average (bulk motion)
    for h = 1:nPhases,
        mnb_P(h) = mean(nonzeros(im_ph_P(:,:,4,h)));
        im_ph_P(:,:,:,h) = im_ph_P(:,:,:,h) - mnb_P(h);
        mnb_M(h) = mean(nonzeros(im_ph_M(:,:,4,h)));
        im_ph_M(:,:,:,h) = im_ph_M(:,:,:,h) - mnb_M(h);
        mnb_S(h) = mean(nonzeros(im_ph_S(:,:,4,h)));
        im_ph_S(:,:,:,h) = im_ph_S(:,:,:,h) - mnb_S(h);
    end

    % subtract off direction-averaged mean phase
    for ph = 1:nPhases,
        im_phase(:,:,:,ph,1) = im_ph_P(:,:,:,ph)-mn_P;
        im_phase(:,:,:,ph,2) = im_ph_M(:,:,:,ph)-mn_M;
        im_phase(:,:,:,ph,3) = im_ph_S(:,:,:,ph)-mn_S;
    end
    for ph = 1:nPhases,
        for sl = 1:nSlices, 
            im_mag(:,:,sl,ph) = mean(im_mag(:,:,sl,ph,:),5);
        end
    end
    
    im_unwrap = im1a_ph;
%------------------------------------4D-----------------------------------%
else
  switch method
    case 'Slow'
      w=waitbar(0,'Progress');
      tic
      for slice = 1:nSlices,
          for phase = 1:nPhases,
              im_unwrap(:,:,slice,phase) = QualityGuidedUnwrap2D_r1(...
                                squeeze(im_mag(:,:,slice,phase)),...
                                squeeze(im_phase(:,:,slice,phase)));
  %             im_unwrap(:,:,slice,phase) = GoldsteinUnwrap2D_r1(...
  %                               squeeze(im_mag(:,:,slice,phase)),...
  %                               squeeze(im_phase(:,:,slice,phase)));
              stat=sprintf('Slice: %g/%g, Phase: %g/%g',slice,nSlices,phase,nPhases);
              waitbar((slice-1+(phase/nPhases))/nSlices,w,stat)
          end
      end
      delete(w)
      runtime=toc;
    case 'Fast'
      tic
      im_unwrap = fastUnwrap(im_mag,im_phase);
      runtime=toc;
  end
    mean_phase = mean(im_unwrap,4);
    for k = 1:nPhases,
        im_unwrap(:,:,:,k) = im_unwrap(:,:,:,k) - mean_phase(:,:,:);
    end
end
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


function bulkmotionBtn_Callback(hObject, eventdata, handles)
% Call bulk motion calculation function & report result

% Conversion factors
% microns per radian of phase
scale = str2num(get(handles.conversionScaleBox,'String'));

% Set up
P = handles.images{handles.curImg}.phase;

% Call calculation function
direction = handles.direction;
showplots = true;
amplitude = bulkmotion(P,scale,direction,showplots);

% Report back to user
message = sprintf('Bulk motion amplitude: %.2g microns',amplitude);
statusMsg(handles,message)
guidata(hObject,handles);


function conversionScaleBox_Callback(hObject, eventdata, handles)


function conversionScaleBox_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function helpBtn_ClickedCallback(hObject, eventdata, handles)
msgbox(['Contact Mark Wagshul <mark.wagshul@einstein.yu.edu> ',...
  'or Alex Krolick <amk283@cornell.edu>'],'Help')
help('MREView')
