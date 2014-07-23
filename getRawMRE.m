function [mag, phase] = getRawMRE()

% Extract data from imaging formats such as DICOM and NIFTI into MATLAB
% arrays.
%   
% Outputs:
%   im_mag:    4- or 5-D Matrix (x,y,z-slice,phase,[direction])
%   im_phase:  4- or 5-D Matrix (x,y,z-slice,phase,[direction])
%
% Dependencies:
%   (+) NIFTI Toolbox: /gmrrc/mrbin/GMRRC/NIFTI
%   (+) BET Toolbox: /gmrrc/mrbin/fsl.cver/bin/bet
%   (*) dcmdump: /gmrrc/mrbin/dcmdump
%   (*) GNU core utils
%    *  Optional, but faster than using MATLAB's DICOM toolbox for metadata
%    +  Required for 5D DICOM processing
%
% Authors:
%   Mark Wagshul <mark.wagshul@einstein.yu.edu>
%   Alex Krolick <amk283@cornell.edu>
%
% See also getMRESinkus, MREAnalyser

[f p index] = uigetfile({...
    '*.dicom','DICOM MRE image';...
    '*.dicom','DICOM Motion-encoded MRE image';...
    '*.nii',  'NIFTI image'});

switch index
  case 1 % MRE Image
    [mag,phase] = readDICOM4D(p,f);
  case 2 % Motion-sensitized MRE image
    [mag,phase] = readDICOM5D(p,f);
  case 3 % NIFTI Image
    [mag,phase] = readNIFTI(p,f);
end

function [mag,phase] = readDICOM4D(p,f)
  info = dicominfo([p f]);
  nSlices = info.Private_2001_1018;
  nPhases = info.Private_2001_1081;
  % Alternately, use dcmdump (only works on MRRC cluster):
  % [~,result] = system(['dcmdump ' p f ' | grep "(2001,1018)" | awk ''{print $3}''']);
  % nSlices = str2num(result);
  % [~,result] = system(['dcmdump ' p f ' | grep NumberOfTemporalPositions | awk ''{print $3}'' | head -n 1 | sed ''s/\[//g'' | sed ''s/]//g''']);
  % nPhases = str2num(result);
  im = double(squeeze(dicomread([p f])));
  for j = 1:2*nSlices,
    for k = 1:nPhases,
      im_(:,:,j,k) = im(:,:,(j-1)*nPhases+k);
    end
  end

  mag = im_(:,:,1:nSlices,:);
  phase = im_(:,:,nSlices+1:2*nSlices,:);
  phase = double(phase)*pi/2048-pi;

  
function [mag,phase] = readDICOM5D(p,f)
  filename = [p '../RAW/' strtok(f,'.') '.nii'];
  if exist(filename) > 0,
    im = lunii('Select nifti image',filename);
  else
    filename = [p '../RAW/' strtok(f,'.') '.nii.gz'];
    im = lunii('Select nifti image',filename);
  end
  im = im.img;
  brainFilename = [p '../RAW/' strtok(f,'.') 'Brain.nii.gz'];
  if exist(brainFilename) > 0,
    imBrain = lunii('Select NIFTI BET image',brainFilename);
  else
    imBrain = lunii('Select NIFTI BET image','');
  end
  imBrain = imBrain.img;
  
  info = dicominfo([p f]);
  nSlices = info.Private_2001_1018;
  nDirs   = info.Private_2001_1081;
  % [~,result] = system(['dcmdump ' p f ' | grep NumberOfTemporalPositions | awk ''{print $3}'' | head -n 1 | sed ''s/\[//g'' | sed ''s/]//g''']);
  % nDirs = str2num(result);
  nPhases = size(im,4) / 2 / nDirs;
  k1=1;
  for dir = 1:nSlices, %nDirs
    for ph = 1:nPhases,
      mag(:,:,:,ph,dir) = im(:,:,:,k1).* int16(imBrain>0);
      phase(:,:,:,ph,dir) = im(:,:,:,size(im,4) / 2 + k1).* int16(imBrain>0);
      k1 = k1 + 1;
     end
  end
  phase = double(phase)*pi/2048-pi;

  
function [mag, phase] = readNIFTI(p,f)
  im = load_untouch_nii([p f]);
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
