function [im_mag, im_phase] = getMREimages(sliceRange,showFig,P_axes)

% Extract magnitude and phase images from DICOM files.
%
% Inputs:
%   sliceRange:  Range, e.g., 1:8
%   showFig:     Display images, true or false (faster with false)
%
% Outputs:
%   im_mag:      4-D Matrix (x,y,z-slice,phase)
%   im_phase:    4-D Matrix (x,y,z-slice,phase)
%
% Depends:
%   dcmdump: /gmrrc/mrbin/dcmdump
%   GNU core utils
%
% Authors:
%   Mark Wagshul <mark.wagshul@einstein.yu.edu>
%   Alex Krolick <amk283@cornell.edu>
%
% See also getMRESinkus

[f p] = uigetfile('*.dicom');
im1 = double(squeeze(dicomread([p f])));

[status,result] = system(['dcmdump ' p f ' | grep "(2001,1018)" | awk ''{print $3}''']);
nSlices = str2num(result);
[status,result] = system(['dcmdump ' p f ' | grep NumberOfTemporalPositions | awk ''{print $3}'' | head -n 1 | sed ''s/\[//g'' | sed ''s/]//g''']);
nPhases = str2num(result);

for j = 1:2*nSlices,
    for k = 1:nPhases,
        im1a(:,:,j,k) = im1(:,:,(j-1)*nPhases+k);
    end
end

im1a_r = im1a(:,:,1:nSlices,:);
im1a_i = im1a(:,:,nSlices+1:2*nSlices,:);

centPos = [];
%for k = 1:nSlices, 
for j = sliceRange,
    [im1a_ph(:,:,:,j),centPos] = nnUnwrap(squeeze(im1a_r(:,:,j,:)),squeeze(im1a_i(:,:,j,:)),1,centPos,showFig,P_axes);
end

%for j = 1:nSlices,
for j = sliceRange,
    for k = 1:nPhases,
        im_phase(:,:,j,k) = im1a_ph(:,:,k,j);
    end
end

im_mag = im1a_r;

mean_phase = mean(im_phase,4);
for k = 1:nPhases,
    im_phase(:,:,:,k) = im_phase(:,:,:,k) - mean_phase(:,:,:);
end
