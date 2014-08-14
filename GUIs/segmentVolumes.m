function varargout = segmentVolumes(varargin)
% Isolate and calculate volume of structures in MRI images.
% 
% Opening files 
% -------------
% Open a NIFTI or DICOM image by clicking "Open". MRI2MAT returns the image
% data and metadata in Matlab matrix format. Coronal, Axial, and Sagittal
% slices are shown in the XZ, XY, and YZ preview planes, respectively. (If
% not, try preprocessing the images with `dcm2nii` and FSL's `bet` and
% `fast` programs or editing the matrices directly in Matlab.) A 3D preview
% of the brain sliced along the sagittal plane is also shown in a separate
% window, which may be closed without causing problems. Note: If the file
% is 4D, the 'Series' slider changes which volume (3D) is shown. Operations 
% like `Isolate` apply only to the active Series.
% 
% Navigation
% ----------
% Left-clicking in any image plane adjusts the other image planes so that
% they are aligned with that coordinate. E.g., clicking (y,z) = (10,20) in
% the sagittal plane shows slice 20 in the axial window and a slice with
% y=10 in the coronal plane. If open file has more than three dimensions,
% for example if the NIFTI file contains a breakdown of CSF, gray matter,
% and white matter as separate volumes, the "Series" slider changes
% th active volume (affects 4D images only).
% 
% Isolating Regions by Connectivity
% ---------------------------------
% The volume can be reduced to include only specific regions by by clicking
% "Isolate." The program segments the selected volume by connectivity,
% finding continous non-zero elements in the matrix and assigning labels to
% each discrete body. The user is shown a view of the axial plane with each
% of these regions color-coded. Clicking a region adds it to the selection.
% Selection continues until the user presses Enter.
% 
% How Volume is Calculated
% ------------------------
% Every time the view updates, volume is recalculated. The NIFTI image
% header is read to determine the units (m, cm, mm) and dimensions of the
% voxels. The volume is calculated as the sum of all non-zero voxels in the
% matrix.
% 
% Select/Cut
% ----------
% To separate the desired region from other connected bodies, you can
% manually remove pieces of the image. To do so, right-click 4 points which
% define a box. The first 2 points determine the Y and Z extents of the
% cut; the second 2 determine the extents in X. Think of this as describing
% a rectangle in the sagittal plane and then setting its width in the
% coronal or axial planes. The action taken depends on the value of
% "Selection":
% 
%   "Cut"    -- Everything inside the box is set to 0
%   "Select" -- Everything outside the box is set to 0
% 
% Erode/Dilate
% ------------
% The "Erode" and "Dilate" buttons execute the Matlab Image Processing
% Toolbox functions `imerode` and `imdilate`. Erosion can be useful to
% break small connections between regions, but should be used sparingly
% because the effect can be quite destructive. Dilation is the inverse
% operation of erode, but doesn't necessarily 'undo' an erosion.%
% 
% Undo
% ----
% Pressing the "Backspace" key clears the last point from the cut
% selection. The "<<" button clears the entire active cut selection. Note:
% There isn't a function to reset a cut that has been already made.
% 
% 3D View
% -------
% The "3D" button pops open a new window containing a 3D isosurface of the
% current volume.
% 
% Saving
% ------
% Clicking the "Save" button in the toolbar saves the image matrix in a
% .mat-file along with the image header (metadata) and the volume with
% units. The file has the same name and is saved in the same folder as the
% opened image. The path, filename, and volume are added to the .csv file
% in the "Log" text box at the top of the GUI. A warning dialog appears if
% for some reason either file cannot be saved (i.e., a permissions error).
% 
% Dependencies
% ------------
% MRI2MAT: Convert NIFTI and DICOM images to matrix format
% GUIDE:   GUI created using GUIDE toolbox
% Matlab Image Processing Toolbox
% 
% Author:
% Alex Krolick <amk283@cornell.edu>
% 
% See also: MRI2MAT, IMERODE, IMDILATE, segmentVolumes.fig
% 
% Edit the above text to modify the response to help segmentVolumes
% 
% Last Modified by GUIDE v2.5 13-Aug-2014 10:34:05
%
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

% SETUP
handles.M = [];
pt = [1,1,1,1];
handles.pt = pt;
handles.volume = 1;
handles.cut = {[],[],[],[]};
handles.ptSelect =1;

% Parameters for erosion/dilation
rad = 1;
[xx,yy,zz] = ndgrid(-rad:rad);
nhood = sqrt(xx.^2 + yy.^2 + zz.^2) <= rad;
handles.nhood = nhood;


% END SETUP

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
% Open a 3 or 4D file using mri2mat. It should contain at least one
% volume matrix with [x,y,z] aligned in the Right/Anterior/Superior 
% coordinate system. Additional volumes add a dimension: [x,y,z,volume_num]
% In the GUI, the planes are: Axial XY, Coronal XZ, Sagittal YZ.

% Open file
[M,header] = mri2mat();

% Set volume selection slider properties based on number of volumes
% detected
dims = size(M);
if dims(1)>=512 % most likely an up-sampled image
  M = impyramid(M,'reduce'); % downsample by a factor of 2
  % scale pixel dimensions accordingly (see nifti_units_lookup function):
  header.dime.pixdim(2:3) = header.dime.pixdim(2:3)*2; 
  dims = size(M); % get new dimensions
end
if length(dims)==2, dims(3)=1; end % pad dims if only 2D
if length(dims)==3, dims(4)=1; end % pad dims if only 3D
set(handles.volumeSelect,'Max',dims(4)+.01);
set(handles.volumeSelect,'Min',1)
set(handles.volumeSelect,'SliderStep',[1 1]/dims(4));
% `handles.pt` is used to determine slice planes.
% E.g., if `pt` is [65,65,50,2], the sagittal slice is M(65,:,:,2)
% The next line sets the starting point to the middle of the brain.
pt = ceil(dims/2); 
set(handles.volumeSelect,'Value',pt(4));
handles.M = M;
handles.pt = pt;
handles.header = header;
statusMsg(handles,sprintf('%s %s',...
  'Right-click twice in the sagittal plane to outline a cut,',...
  'then twice in the coronal plane to set the width.'))

updateVolume(handles)

guidata(hObject, handles);
update_Axes(hObject,handles);

% 3D image
% Create a floating figure window
handles.floatingFig = figure('name','3D View');
colormap gray
D = M(round(dims(1)/2):end,:,:,1); % split along sagittal plane
Ds = smooth3(D);
hiso = patch(isosurface(Ds,0),...
	'FaceColor',[.5,.5,.8],...
	'EdgeColor','none');
	isonormals(Ds,hiso)
hcap = patch(isocaps(D,0),...
	'FaceColor','interp',...
	'EdgeColor','none');
view(35,30) 
axis tight 
pixdim = header.dime.pixdim(2:4); % get voxel dimensions
daspect(1./[pixdim(1),pixdim(2),pixdim(3)]) % scale proportional to voxel shape
lightangle(90,0);
set(gcf,'Renderer','opengl'); lighting phong
set(hcap,'AmbientStrength',1)
set(hiso,'SpecularColorReflectance',0,'SpecularExponent',50)


function isolateButton_Callback(hObject, eventdata, handles)
% Ask user to select which subvolume they want to calculate, then calculate
% volume and display a 3D image of the region. Save the matrices and figure
% to disk.

% Region selection and isolation
header = handles.header;
M = handles.M;
v = handles.pt(4);
M = M(:,:,:,v);

% Label connected components
  L = labelmatrix(bwconncomp(M)); % find and label connected regions
  %L(L(:)==0) = NaN; % set all background (0) values to NaN to avoid counting

  % Sort labeled regions by volume
  s = regionprops(L, {'Area', 'PixelIdxList'});
  areas = cat(1, s.Area);
  [~, sort_order] = sort(areas,'descend');
  s2 = s(sort_order);
  for k = 1:numel(s2)
     pixelids = s2(k).PixelIdxList;
     L(pixelids) = k;
  end

  Lslice = L(:,:,handles.pt(3)); % the active axial slice
  % Ask for selection
  % Show the user `Lslice` and ask them to click on the colored region they
  % want to use for the volume calculation
  figure
  colormap jet;
  imagesc(Lslice,[0 15]);
  title 'Select regions to keep. Press enter when finished.'
  coords = round(ginput()); % get label of the region the user wants to see
  temp = zeros(size(M));
  for i = 1:size(coords,1)
    val = Lslice(coords(i,2),coords(i,1)); % get value of the label
    temp(L==val) = M(L==val);
  end
  
M = temp;

statusMsg(handles,'Segment isolation complete')
handles.M(:,:,:,v) = M;
guidata(hObject,handles)
update_Axes(hObject,handles)


function [space_units,time_units] = nifti_units_lookup(code)
% See the NIFTI specification: 
% http://nifti.nimh.nih.gov/pub/dist/src/niftilib/nifti1.h
% Bits 0-2 are spatial, 3-5 temporal, 6-7 unused
code = dec2bin(code);
space_code = bin2dec(code(end-2:end-0));
time_code  = bin2dec(code(:,end-3));
switch space_code
  case 1, space_units = 1E0;  %meters
  case 2, space_units = 1E-3; %millimeters
  case 3, space_units = 1E-6; %micrometers
end

switch time_code
  case 1, time_units = 1E0;  %seconds
  case 2, time_units = 1E-3; %milliseconds
  case 3, time_units = 1E-6; %microseconds
end

% Callbacks for clicks in the view planes (axes1,axes2,axes3):

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
% Sagittal (YZ) Plane
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
% Redraw the view planes. Called after opening a file or when requested by
% a callback function
f = [handles.axes1,handles.axes2,handles.axes3];
M = handles.M;
pt = handles.pt;
colormap gray
axis tight
% Cut out the selected region if selection is complete
% If selection mode is 'Select', cut out everything else
c = handles.cut;
if ~isempty(c{1}) && ~isempty(c{2}) && ~isempty(c{3}) && ~isempty(c{4})
  switch get(handles.selectModeMenu,'Value')
    case 2
      M(min([c{3}(1),c{4}(1)]) : max([c{3}(1),c{4}(1)]),... % X extent
        min([c{1}(2),c{2}(2)]) : max([c{1}(2),c{2}(2)]),... % Y extent
        min([c{1}(3),c{2}(3)]) : max([c{1}(3),c{2}(3)]),... % Z extent
        c{1}(4)) = 0;
      handles.cut = {[],[],[],[]}; %reset selection
    case 1
      Z = zeros(size(M));
      Z(min([c{3}(1),c{4}(1)]) : max([c{3}(1),c{4}(1)]),... % X extent
        min([c{1}(2),c{2}(2)]) : max([c{1}(2),c{2}(2)]),... % Y extent
        min([c{1}(3),c{2}(3)]) : max([c{1}(3),c{2}(3)]),... % Z extent
        c{1}(4)) = 1;
      handles.cut = {[],[],[],[]}; %reset selection
      M(Z==0) = 0;
  end
  handles.M = M; %save changes
  guidata(hObject,handles)
end

% Redraw the images in each plane. Transposes ensure Z points vertically
% rather than horizontally.'HitTest','off' makes clicks pass through to the
% parent axes so mouse click coordinates can be recorded with axes callbacks
imagesc(transpose(squeeze(M(:,pt(2),:,pt(4)))),'Parent',f(1),'HitTest','off');
imagesc(squeeze(M(:,:,pt(3),pt(4))),'Parent',f(2),'HitTest','off');
imagesc(transpose(squeeze(M(pt(1),:,:,pt(4)))),'Parent',f(3),'HitTest','off');
for i=1:3, axis(f(i),'tight'); end
updateCoords(handles)
updateVolume(handles)
guidata(hObject,handles)


function volumeSelect_Callback(hObject, eventdata, handles)
% volume selection slider callback
v = round(get(hObject,'Value')); %get slider position
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
  sprintf('(%g,%g,%g)',x,y,z));
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


function view3DBtn_Callback(hObject, eventdata, handles)
% 3D image of selected region
header = handles.header;
M = handles.M;
v = handles.pt(4);
M = M(:,:,:,v)>0;
figure('name','3D View')
Ds = smooth3(M);
isosurface(Ds,0)
pixdim = header.dime.pixdim(2:4); % get voxel dimensions
daspect(1./[pixdim(1),pixdim(2),pixdim(3)]) % scale proportional to voxel shape

function selectModeMenu_Callback(hObject, eventdata, handles)


function selectModeMenu_CreateFcn(hObject, eventdata, handles)

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function erodeBtn_Callback(hObject, eventdata, handles)
M = handles.M;
nhood = handles.nhood;
v = handles.pt(4);
M = M(:,:,:,v);
M = imerode(M,nhood);
statusMsg(handles,'Erosion complete')
handles.M(:,:,:,v) = M;
guidata(hObject,handles)
update_Axes(hObject,handles)

function dilateBtn_Callback(hObject, eventdata, handles)
M = handles.M;
nhood = handles.nhood;
v = handles.pt(4);
M = M(:,:,:,v);
M = imdilate(M,nhood);
statusMsg(handles,'Dilation complete')
handles.M(:,:,:,v) = M;
guidata(hObject,handles)
update_Axes(hObject,handles)

function volumeTxt_ButtonDownFcn(hObject, eventdata, handles)
updateVolume(handles);

function updateVolume(handles)
% calculate volume of active region

header = handles.header;
M = handles.M;
v = handles.pt(4);
M = M(:,:,:,v)>0;

% Volume calculation
[space_units,~] = nifti_units_lookup(header.dime.xyzt_units);
pixdim = handles.header.dime.pixdim(2:4); % voxel dimensions in real units
voxels = sum(M(:)>0); % add up total voxels in selected region
volume = voxels*pixdim(1)*pixdim(2)*pixdim(3)*space_units^3; % vox -> cubic meters

% Display volume
volume_mL = volume*1E6; % cubic meters -> milliters
volumeStr = sprintf('%.1f mL',volume_mL); 
%fprintf([volumeStr '\n']);

% Update volume text box
set(handles.volumeTxt,'String',volumeStr);
pause(.01)


function saveBtn_ClickedCallback(hObject, eventdata, handles)
% Save GUI data to disk

% Inputs
header = handles.header;
M = handles.M;
volume = get(handles.volumeTxt,'String'); % e.g., '0.00 mL'

% Outputs
% - save a .mat file next to the original data file
% - append metadata and volume to a log file in the current directory
matfile = [header.path strtok(header.filename,'.') '.mat']; 
csvname = get(handles.csvFileBox,'String');
csvfile = fopen(csvname,'a'); % open with append-only permissions

% Write files
try
  fprintf(csvfile,['\n' '"' header.path '","' header.filename '","' volume '"']);
catch
  errordlg(['Couldn''t write to ' csvname])
  errflag = true;
end
try
  save(matfile, 'header', 'M', 'volume');
catch
   errordlg(['Couldn''t write to ' matfile])
   errflag = true;
end

% Close
fclose(csvfile);

% Update dialog
switch exist('errflag','var')
  case true
    statusMsg(handles,'Save failed.')
  case false
    statusMsg(handles,'Saved.')
end


function csvFileBox_Callback(hObject, eventdata, handles)


function csvFileBox_CreateFcn(hObject, eventdata, handles)
set(hObject,'String',[pwd '/Volumes.csv']);
