function varargout = MREAnalyser(varargin)
% MREANALYSER Preview and phase unwrap MRE sequences.
%
% See also: GUIDE, GUIDATA, GUIHANDLES, UNWRAPPER, MRI2MAT
% 
% Authors:
% Alex Krolick <amk283@cornell.edu>
% Mark Wagshul <mark.wagshul@einstein.yu.edu>

% Edit the above text to modify the response to help MREAnalyser

% Last Modified by GUIDE v2.5 16-Jul-2014 12:23:14

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
statusMsg(handles,'Open an image file')

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

[f p index] = uigetfile({...
    '*.dicom','DICOM MRE image';...
    '*.dicom','DICOM Motion-encoded MRE image';...
    '*.nii*',  'NIFTI image'});

w = handles.sliceWindow;
i = handles.curImg;

% info = dicominfo(f);
% nSlices = double(info.Private_2001_1018);
% nPhases = double(info.Private_2001_1081);

statusMsg(handles,'Opening image...');

switch index
    case 1 % MRE Image
        im = double(squeeze(dicomread([p f])));
        statusMsg(handles,'Reading metadata...');
        [~,result] = system(['dcmdump ' p f ' | grep "(2001,1018)" | awk ''{print $3}''']);
        nSlices = str2num(result);
        [~,result] = system(['dcmdump ' p f ' | grep NumberOfTemporalPositions | awk ''{print $3}'' | head -n 1 | sed ''s/\[//g'' | sed ''s/]//g''']);
        nPhases = str2num(result); 

        for j = 1:2*nSlices,
            for k = 1:nPhases,
                im1a(:,:,j,k) = im(:,:,(j-1)*nPhases+k);
            end
        end
        mag = im1a(:,:,1:nSlices,:);
        phase = im1a(:,:,nSlices+1:2*nSlices,:);
        phase = double(phase - 2048)*pi/2048;
        %[mag, phase] = getMREimages(1:30,false,w,f,p);
        handles.images{i}.mag = mag;
        handles.images{i}.phase = phase;
        handles.activeImg = phase(:,:,:,:,1);
        set(handles.dirBtnGrp,'Visible','off');  
    case 2 % Motion-sensitized MRE image
        %[mag, phase] = getMRESinkus(false,f,p,handles.sliceWindow);
        filename = [p '../RAW/' strtok(f,'.') '.nii'];
        if exist(filename) > 0,
            im = lunii('Select nifti image',filename);
        else
            filename = [p '../RAW/' strtok(f,'.') '.nii.gz'];
            im = lunii('Select nifti image',filename);
        end
        im = im.img;
        nSlices = size(im,3);
        brainFilename = [p '../RAW/' strtok(f,'.') 'Brain.nii.gz'];
        if exist(brainFilename) > 0,
            imBrain = lunii('Select NIFTI BET image',brainFilename);
        else
            imBrain = lunii('Select NIFTI BET image','');
        end
        imBrain = imBrain.img;
        [~,result] = system(['dcmdump ' p f ' | grep NumberOfTemporalPositions | awk ''{print $3}'' | head -n 1 | sed ''s/\[//g'' | sed ''s/]//g''']);
        nDirs = str2num(result);
        nPhases = size(im,4) / 2 / nDirs;
        k1=1;
        for dir = 1:nDirs,
            for ph = 1:nPhases,
                mag(:,:,:,ph,dir) = im(:,:,:,k1).* int16(imBrain>0);
                phase(:,:,:,ph,dir) = im(:,:,:,size(im,4) / 2 + k1).* int16(imBrain>0);
                k1 = k1 + 1;
            end
        end
        phase = double(phase - 2048)*pi/2048;
        handles.images{i}.mag = mag;
        handles.images{i}.phase = phase;
        handles.activeImg = phase(:,:,:,:,1);
        set(handles.dirBtnGrp,'Visible','on');
    case 3 % NIFTI
        im = load_untouch_nii([p f]);
        mag = im.img;
        phase = []; 
        [nX,nY,nSlices,nPhases,nDirs] = size(mag);
        handles.images{i}.mag = mag;
        handles.images{i}.phase = phase;
        handles.activeImg = mag(:,:,:,:);
        set(handles.dirBtnGrp,'Visible','off');
        set(handles.magToggle,'Visible','off');
end

if index==1 || index==2
    statusMsg(handles,'Setting up workspace...');

    p = handles.phaseSelect;
    s = handles.sliceSelect;
    set(s,'Min',1), set(s,'Max',nSlices), set(s, 'Value',1);
    set(p,'Min',1), set(p,'Max',nPhases), set(p, 'Value',1);
    set(s, 'SliderStep', [1 1]/max((nSlices-1),1));
    set(p, 'SliderStep', [1 1]/max((nPhases-1),1));

    guidata(hObject,handles);
    set(handles.magToggle,'SelectedObject',handles.phaseBtn);
    updateCAxis(handles);
    redraw(handles);
    statusMsg(handles,'Image loaded')
else
    statusMsg(handles,'Setting up workspace...');

    p = handles.phaseSelect;
    s = handles.sliceSelect;
    set(s,'Min',1), set(s,'Max',nSlices), set(s, 'Value',1);
    set(p,'Min',1), set(p,'Max',nPhases), set(p, 'Value',1);
    set(s, 'SliderStep', [1 1]/max((nSlices-1),1));
    set(p, 'SliderStep', [1 1]/max((nPhases-1),1));

    set(handles.magToggle,'SelectedObject',handles.magBtn);
    guidata(hObject,handles);
    updateCAxis(handles);
    redraw(handles);
    statusMsg(handles,'Image loaded')
end

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
dims = size(handles.images{i}.phase);
nSlices= dims(3); nPhases= dims(4);
im_mag = handles.images{i}.mag;
im_phase = handles.images{i}.phase;
im_unwrap = im_phase;
%-------------------------5D-------------------------%
if length(dims) == 5, 
%     im1a_ph = im_phase;
%     for dir = 1:dims(5);
%         for slice = 1:nSlices,
%             for phase = 1:nPhases,
%                 im1a_ph(:,:,slice,phase) = QualityGuidedUnwrap2D_r1(...
%                                   squeeze(im_mag(:,:,slice,phase,dir)),...
%                                   squeeze(im_phase(:,:,slice,phase,dir)));
%                 fprintf('Slice: %g/%g, Phase: %g/%g, Dir: %g\n',slice,nSlices,phase,nPhases,dir) 
%             end
%         end
%     end
    im1a_ph = unwrapper(im_mag,im_phase);
     % Subtract phase due to background field 
    im_ph_P = im1a_ph(:,:,:,:,1)-im1a_ph(:,:,:,:,4); % Phase direction
    im_ph_M = im1a_ph(:,:,:,:,2)-im1a_ph(:,:,:,:,4); % Magnitude direction
    im_ph_S = im1a_ph(:,:,:,:,3)-im1a_ph(:,:,:,:,4); % Slice direction

    % Time (phase) average
    mn_P = mean(im_ph_P,4);
    mn_M = mean(im_ph_M,4);
    mn_S = mean(im_ph_S,4);

    % Direction average (bulk motion)
    %size(im_ph_P(4))
    for h = 1:8,
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
%     for slice = 1:nSlices,
%         for phase = 1:nPhases,
%             im_unwrap(:,:,slice,phase) = QualityGuidedUnwrap2D_r1(...
%                               squeeze(im_mag(:,:,slice,phase)),...
%                               squeeze(im_phase(:,:,slice,phase)));
% %             im_unwrap(:,:,slice,phase) = GoldsteinUnwrap2D_r1(...
% %                               squeeze(im_mag(:,:,slice,phase)),...
% %                               squeeze(im_phase(:,:,slice,phase)));
%             fprintf('Slice: %g/%g, Phase: %g/%g\n',slice,nSlices,phase,nPhases) 
%         end
%     end
    im_unwrap = unwrapper(im_mag,im_phase);
    mean_phase = mean(im_unwrap,4);
    for k = 1:nPhases,
        im_unwrap(:,:,:,k) = im_unwrap(:,:,:,k) - mean_phase(:,:,:);
    end
end

handles.images{i}.phase = im_unwrap;
handles.images{i}.mag = im_mag;
handles.activeImg = handles.images{i}.phase(:,:,:,:,1);
updateCAxis(handles);
guidata(hObject,handles);
redraw(handles);
