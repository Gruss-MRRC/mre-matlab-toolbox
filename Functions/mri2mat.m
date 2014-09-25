function [mag, varargout] = mri2mat()
% Convert NIFTI and DICOM images to matrix format.
%
% Inputs:
%   .nii, .nii.gz, .dcm, .dicom files (specified interactively)
%   
% Outputs:
%   mag:     4 or 5-D matrix (x,y,slice,phase,[direction])
%  [phase]:  4 or 5-D matrix (x,y,slice,phase,[direction]) *optional
%  [info]:   structure containing metadata about the file  *optional
%
% Dependencies:
%   + NIFTI Toolbox: /gmrrc/mrbin/GMRRC/NIFTI
%   ~ BET Toolbox: /gmrrc/mrbin/fsl.cver/bin/bet
%   * dcmdump: /gmrrc/mrbin/dcmdump
%   * GNU core utils (grep, awk, etc.)
%   --------------
%   +  Required to work with NIFTI images
%   ~  Required to create masks, but not called in this program
%   *  Optional, but faster than using MATLAB's DICOM toolbox for metadata
% 
% Authors:
%   Mark Wagshul <mark.wagshul@einstein.yu.edu>
%   Alex Krolick <amk283@cornell.edu>
%
% See also getMRESinkus, MRE_Preview, nnUnwrap, load_untouch_nii, load_nii,
% dicomread, dicominfo

% Open file
[f p index] = uigetfile({...
    '*.dicom', '3D DICOM to 4D Mag & Phase Matrices';...
    '*.dicom', '3D DICOM to 5D Mag & Phase Matrices';...
    '*.dicom', '3D DICOM to 5D Mag & Phase Matrices with NIFTI mask';...
    '*.nii',   '4D NIFTI to 5D Mag & Phase Matrices';...
    '*.nii',   'NIFTI File';...
    '*.nii.gz','NIFTI Archive';...
    '*.mat',   'MAT File from MRE_Preview';});
  
% Do something based on filetype
switch index
  case 1 % MRE Image
    [mag,phase,info] = readDICOM4D(p,f);
    no_phase = false;
  case 2 % Motion-sensitized MRE image
    [mag,phase,info] = readDICOM5D(p,f);
    no_phase = false;
  case 3 % Motion-sensitized MRE image + mask
    [mag,phase,info] = readDICOM5D_Mask(p,f);
    no_phase = false;
  case 4 % 5D NIFTI
    [mag,phase,info] = readNIFTI5D(p,f);
    no_phase = false;
  case 5 % NIFTI file
    [mag,info] = readNIFTI(p,f);
    no_phase = true;
  case 6 % NIFTI archive
    [mag,info] = readNIFTI(p,f);
    no_phase = true;
  case 7 % MAT file
    vars = open([p f]);
    mag = vars.M;
    phase = vars.P;
    info = vars.header;
    no_phase = false;
end

% Add data about which file was opened to header info
info.filename = f;
info.path = p;

% Process optional outputs (phase and info)
if nargout == 3, 
  varargout{1} = phase; 
  varargout{2} = info; % return mag, phase, and info
elseif nargout == 2 && no_phase, 
  varargout{1} = info; % return mag and info
elseif nargout == 2 && ~no_phase,
  varargout{1} = phase; % return mag and phase
end


function [mag,phase,info] = readDICOM4D(p,f)
% based on getMREimages.m

  if exist('/gmrrc/mrbin/dcmdump','file')  && isunix()
    % Use `dcmdump` to get metadata (only works on MRRC cluster):
    [~,result] = system(['dcmdump ' p f ...
      ' | grep "(2001,1018)" | awk ''{print $3}''']);
    nSlices = str2num(result);
    [~,result] = system(['dcmdump ' p f...
      ' | grep NumberOfTemporalPositions'...
      ' | awk ''{print $3}'' | head -n 1 | sed ''s/\[//g'' | sed ''s/]//g''']);
    nPhases = str2num(result);
    info = struct();
  else
    % Fall back to using Matlab utilities
    info = dicominfo([p f]);
    nSlices = info.Private_2001_1018;
    nPhases = info.Private_2001_1081;
  end
  
  im = dicomread([p f]);
  im = double(squeeze(im));
  
  for j = 1:2*nSlices,
    for k = 1:nPhases,
      im_(:,:,j,k) = im(:,:,(j-1)*nPhases+k);
    end
  end

  mag = im_(:,:,1:nSlices,:);
  phase = im_(:,:,nSlices+1:2*nSlices,:);
  phase = double(phase)*pi/2048-pi;
  

  function [mag,phase,info] = readDICOM5D(p,f)
% Take a 3D DICOM file and output 5D matrices for phase and magnitude by
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
  im = double(squeeze(im));
  
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
  
  
function [mag,phase,info] = readDICOM5D_Mask(p,f)
% Same as READDICOM5D, but mask the output. Usually the mask is created
% using `Brain Extraction Toolbox`, called using 'bet' at the command line.
% See also getMRESinkus.m, readNIFTI5D, readDICOM5D

  % Acquire NIFTI image
  filename = [p '../RAW/' strtok(f,'.') '.nii'];
  if exist(filename) > 0,
    im = lunii('Select nifti image',filename);
  else
    filename = [p '../RAW/' strtok(f,'.') '.nii.gz'];
    im = lunii('Select nifti image',filename);
  end
  im = im.img;
  
  % Acquire NIFTI mask
  brainFilename = [p '../RAW/' strtok(f,'.') 'Brain.nii.gz'];
  if exist(brainFilename) > 0,
    imBrain = lunii('Select NIFTI BET image',brainFilename);
  else
    imBrain = lunii('Select NIFTI BET image','');
  end
  imBrain = imBrain.img;
  mask = int16(imBrain>0);
  
  % Metadata
  if exist('/gmrrc/mrbin/dcmdump','file') && isunix()
    nSlices = size(im,3);
    [~,result] = system(['dcmdump ' p f ...
      ' | grep NumberOfTemporalPositions'...
      ' | awk ''{print $3}'' | head -n 1 | sed ''s/\[//g'' | sed ''s/]//g''']);
    nDirs = str2num(result);
    info = struct();
  else
    info = dicominfo([p f]);
    %nSlices = info.Private_2001_1018;
    nDirs   = info.Private_2001_1081;
  end
  nVolumes = size(im,4);
  nPhases = nVolumes / 2 / nDirs;
  
  % Rearrange volumes
  k1=1;
  for dir = 1:nDirs,
    for ph = 1:nPhases,
      mag(:,:,:,ph,dir) = im(:,:,:,k1).* mask;
      phase(:,:,:,ph,dir) = im(:,:,:,nVolumes / 2 + k1).* mask;
      k1 = k1 + 1;
     end
  end
  
  % Convert to radians
  phase = double(phase)*pi/2048-pi;
  
function [mag,phase,info] = readNIFTI5D(p,f)
% Take a 4D NIFTI file and output 5D matrices for phase and magnitude by
% shuffling the 4th dimension. Assume that the input format is [x,y,z,t],
% where the time axis is filled according to the sequence
% [(mag,dir1,t1),(mag,dir1,t2)...(phase,dirN,tN)]. Decompose into a matrix
% for phase and a matrix for magnitude, which should each be of the form
% [x,y,z,phase,direction], where phase is a temporal dimension. The
% difference between READNIFTI5D and READDICOM5D is that the DICOM header
% contains information about the number of dimensions, but the NIFTI header
% does not. The number of directions is assumed to be 4 in READNIFTI5D.

  im = load_untouch_nii([p f]);
  info = im.hdr;
  im = im.img;
  nVolumes = info.dime.dim(5);
  nDirs=4; % x y z b0
  nPhases = nVolumes / 2 / nDirs;
  k1=1;
  for dir = 1:nDirs,
    for ph = 1:nPhases,
      mag(:,:,:,ph,dir) = im(:,:,:,k1);
      phase(:,:,:,ph,dir) = im(:,:,:,nVolumes / 2 + k1);
      k1 = k1 + 1;
     end
  end
  
  % Convert to radians
  phase = double(phase)*pi/2048-pi;

  
function [mag, info] = readNIFTI(p,f)
% Load a NIFTI file (.nii) or archive (.nii.gz).
% Depends on the NIFTI toolbox for MATLAB.
% See also [1],[2] in the index comments at the end of the file

  im = load_untouch_nii([p f]);
  info = im.hdr;
  mag = im.img;
  phase = [];


function [nii_file, fileLocation] = lunii(openText,fileLocation)
% Load NIFTI (.nii) file (untouch)
% See also load_untouch_nii
  if isempty(fileLocation),
      [f p] = uigetfile('*.nii',openText);
      fileLocation = [p f];
  end
  nii_file = load_untouch_nii(fileLocation);
  
  
%----------------------------------INDEX----------------------------------%
%
% [1] http://brainder.org/2012/09/23/the-nifti-file-format
%
% [2] Information on NIFTI dimension field from [1]:
% The field short dim[8] contains the size of the image array. The first
% element (dim[0]) contains the number of dimensions (1-7). If dim[0] is
% not in this interval, the data is assumed to have opposite endianness and
% so, should be byte-swapped (the nifti standard does not specify a
% specific field for endianness, but encourages the use of dim[0] for this
% purpose. The dimensions 1, 2 and 3 are assumed to refer to space (x, y,
% z), the 4th dimension is assumed to refer to time, and the remaining
% dimensions, 5, 6 and 7, can be anything else. The value dim[i] is a
% positive integer representing the length of the i-th dimension.

  
