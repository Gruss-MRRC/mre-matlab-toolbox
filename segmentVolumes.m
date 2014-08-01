function varargout = segmentVolumes(varargin)
% SEGMENTVOLUMES Isolate and calculate volume of structures in MRI images.
% 
% Open an interactive GUI to place cuts and perform segmentation of volumes
% in an MRI sequence. The intended use is to calculate the volume of the
% ventricles in the human brain using T1 MRI images.
%
% Dependencies:
% MRI2MAT: Convert NIFTI and DICOM images to matrix format
% GUIDE:   GUI created using GUIDE toolbox
%
% Author:
% Alex Krolick <amk283@cornell.edu>
%
% See also: GUIDE, GUIDATA, GUIHANDLES, MRI2MAT

% Edit the above text to modify the response to help segmentVolumes

% Last Modified by GUIDE v2.5 30-Jul-2014 14:53:24

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
% uiwait(handles.segmentVolumes);


function varargout = segmentVolumes_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function openButton_Callback(hObject, eventdata, handles)
[M,~,header] = mri2mat();
handles.M = M;

% Set volume selection slider properties based on number of volumes
% detected
dims = size(M);
if length(dims)==3, dims(4)=1; end % pad dims if only 3D
set(handles.volumeSelect,'Max',dims(4)+.01);
set(handles.volumeSelect,'Min',1)
set(handles.volumeSelect,'SliderStep',[1 1]/dims(4));
pt = ceil(dims/2);
set(handles.volumeSelect,'Value',pt(4));
handles.pt = pt;
handles.header = header;
statusMsg(handles,sprintf('%s %s',...
  'Right-click twice in the sagittal plane to outline a cut,',...
  'then twice in the coronal plane to set the width.'))
guidata(hObject, handles);
update_Axes(hObject,handles);

% 3D image
figure
D = M(round(dims(1)/2):end,:,:,1); % split along sagittal plane
Ds = smooth3(D);
hiso = patch(isosurface(Ds,5),...
	'FaceColor',[.5,.5,.8],...
	'EdgeColor','none');
	isonormals(Ds,hiso)
hcap = patch(isocaps(D,5),...
	'FaceColor','interp',...
	'EdgeColor','none');
view(35,30) 
axis tight 
daspect([1,1,1])
lightangle(90,0);
set(gcf,'Renderer','opengl'); lighting phong
set(hcap,'AmbientStrength',1)
set(hiso,'SpecularColorReflectance',0,'SpecularExponent',50)


function calculateButton_Callback(hObject, eventdata, handles)
% Ask user to select which subvolume they want to calculate, then calculate
% volume and display a 3D image of the region

% Region selection
M = handles.M;
v = get(handles.volumeSelect,'Value');
L = labelmatrix(bwconncomp(M(:,:,:,v)));
L(L(:)==0) = NaN; % set all background (0) values to NaN to avoid counting
% disp('<Volume Histogram>')
% group_histogram = histc(L(:),0:max(L(:)))
% [~,index] = max(group_histogram);
% disp('</Volume Histogram>')
Lslice = L(:,:,handles.pt(3));
h = figure; 
colormap jet
imagesc(Lslice,[0 10]) % TODO don't hardcode the color axis
coords = round(ginput(h));
val = Lslice(coords(2),coords(1));

% Volume calculation
[space_units,~] = nifti_units_lookup(header.dime.xyzt_units);
pixdim = header.dime.pixdim(2:4); % voxel dimensions in real units
voxels = sum(L(:)==val); 
volume = voxels*pixdim(1)*pixdim(2)*pixdim(3)*space_units^3; % cubic meters

% Display volume
volumeMsg = sprintf('Volume of selected region: %g mL',volume*1E6); % m^3 -> mL
statusMsg(handles,volumeMsg);
fprintf([volumeMsg '\n']);

% 3D image of selected region
Z = L==val;
figure(h)
D = M(:,:,:,1).*Z;
Ds = smooth3(D);
hiso = patch(isosurface(Ds,5),...
	'FaceColor',[.5,.5,.8],...
	'EdgeColor','none');
	isonormals(Ds,hiso)
view(35,30) 
axis tight 
daspect([1,1,1])
lightangle(90,0);
set(gcf,'Renderer','opengl'); lighting phong
set(hiso,'SpecularColorReflectance',0,'SpecularExponent',50)

function [space_units,time_units] = nifti_units_lookup(code)
% See the NIFTI specification: 
% http://nifti.nimh.nih.gov/pub/dist/src/niftilib/nifti1.h
% Bits 0-2 are spatial, 3-6 temporal, 7-8 unused
code = dec2bin(code);
space_code = bin2dec(code(end-2:end-0));
time_code  = bin2dec(code(end-5,end-3));
switch space_code
  case '1', space_units = 1E0;  %meters
  case '2', space_units = 1E-3; %millimeters
  case '3', space_units = 1E-6; %micrometers
end

switch time_code
  case '1', time_units = 1E0;  %seconds
  case '2', time_units = 1E-3; %milliseconds
  case '3', time_units = 1E-6; %microseconds
end



function axes1_ButtonDownFcn(hObject, eventdata, handles)
% Coronal (XZ) plane
newpt = get(hObject,'CurrentPoint');
p = handles.segmentVolumes;
oldpt = handles.pt;
switch get(p,'SelectionType')
  case 'normal'%left mouse button click
    handles.pt = round([newpt(1),oldpt(2),newpt(3),oldpt(4)]);   
  case 'alt'%right mouse button click
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


function axes2_ButtonDownFcn(hObject, eventdata, handles)
% Axial (XY) Plane
newpt = get(hObject,'CurrentPoint');
p = handles.segmentVolumes;
oldpt = handles.pt;
switch get(p,'SelectionType')
  case 'normal'%left mouse button click
    
    handles.pt = round([newpt(3),newpt(1),oldpt(3),oldpt(4)]);   
  case 'alt'%right mouse button click
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


function axes3_ButtonDownFcn(hObject, eventdata, handles)
% Sagittal (YZ) Plane, transposed so Z is up
newpt = get(hObject,'CurrentPoint');
p = handles.segmentVolumes;
oldpt = handles.pt;
switch get(p,'SelectionType')
  case 'normal'%left mouse button click
        handles.pt = round([oldpt(1),newpt(1),newpt(3),oldpt(4)]);   
  case 'alt'%right mouse button click
    i = handles.ptSelect;
    handles.cut{i} = round([oldpt(1),newpt(1),newpt(3),oldpt(4)]);
    if handles.ptSelect<5
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

% Cut out the selected region if selection is complete
c = handles.cut;
if ~isempty(c{1}) && ~isempty(c{2}) && ~isempty(c{3}) && ~isempty(c{4})
  M(min([c{3}(1),c{4}(1)]) : max([c{3}(1),c{4}(1)]),... % X extent
    min([c{1}(2),c{2}(2)]) : max([c{1}(2),c{2}(2)]),... % Y extent
    min([c{1}(3),c{2}(3)]) : max([c{1}(3),c{2}(3)]),... % Z extent
    c{1}(4)) = 0;
  handles.cut = {[],[],[],[]};
  handles.M = M;
  guidata(hObject,handles)
end
imagesc(transpose(squeeze(M(:,pt(2),:,pt(4)))),'Parent',f(1),'HitTest','off');
imagesc(squeeze(M(:,:,pt(3),pt(4))),'Parent',f(2),'HitTest','off');
imagesc(transpose(squeeze(M(pt(1),:,:,pt(4)))),'Parent',f(3),'HitTest','off');
for i=1:3, axis(f(i),'tight'); end
updateCoords(handles)
guidata(hObject,handles)


function volumeSelect_Callback(hObject, eventdata, handles)
v = round(get(hObject,'Value'));
handles.pt(4) = v;
statusMsg(handles,sprintf('Changed to volume %g',v))
guidata(hObject,handles)
update_Axes(hObject,handles)


function volumeSelect_CreateFcn(hObject, eventdata, handles)
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function statusMsg(handles,msg)
% Sets status bar text
set(handles.statusBar,'String',msg);
pause(.01)


function updateCoords(handles)
% Displays coordinate values based on the current view
h= handles;
x = h.pt(1);
y = h.pt(2);
z = h.pt(3);
set(h.coords,'String',...
  sprintf('%g,%g,%g',x,y,z));
pause(.01)


function undoButton_Callback(hObject, eventdata, handles)
% Reset the cut selection and selection counter
handles.cut = {[],[],[],[]};
statusMsg(handles,'Selection cleared.')
handles.ptSelect=1;
guidata(hObject,handles)


function segmentVolumes_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to segmentVolumes (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
switch eventdata.Key
  case 'backspace'
    p = handles.ptSelect;
    handles.cut{p} = [];
    handles.ptSelect = handles.ptSelect-1;
    statusMsg(handles,'Undid click.')
end
guidata(hObject,handles)
