function [im_mag, im_phase] = getMRESinkus(showFig)

[f p] = uigetfile('*.dicom','Select P image');
im_p = double(squeeze(dicomread([p f])));

[f p] = uigetfile([p '/*.dicom'],'Select M image');
im_m = double(squeeze(dicomread([p f])));

[f p] = uigetfile([p '/*.dicom'],'Select S image');
im_s = double(squeeze(dicomread([p f])));

[f p] = uigetfile([p '/*.dicom'],'Select Z image');
im_z = double(squeeze(dicomread([p f])));

[status,result] = system(['dcmdump ' p f ' | grep "(2001,1018)" | awk ''{print $3}''']);
nSlices = str2num(result);
[status,result] = system(['dcmdump ' p f ' | grep NumberOfTemporalPositions | awk ''{print $3}'' | head -n 1 | sed ''s/\[//g'' | sed ''s/]//g''']);
nPhases = str2num(result);

nDirs = 4;

k1=1;
for sl = 1:nSlices,
    for ph = 1:nPhases,
        im1a_r(:,:,sl,ph,1)=im_p(:,:,k1);
        im1a_i(:,:,sl,ph,1)=im_p(:,:,size(im_p,3) / 2 + k1);
        im1a_r(:,:,sl,ph,2)=im_m(:,:,k1);
        im1a_i(:,:,sl,ph,2)=im_m(:,:,size(im_p,3) / 2 + k1);
        im1a_r(:,:,sl,ph,3)=im_s(:,:,k1);
        im1a_i(:,:,sl,ph,3)=im_s(:,:,size(im_p,3) / 2 + k1);
        im1a_r(:,:,sl,ph,4)=im_z(:,:,k1);
        im1a_i(:,:,sl,ph,4)=im_z(:,:,size(im_p,3) / 2 + k1);
        k1 = k1 + 1;
    end
end

centPos = [];
for k = 1:nDirs, 
    for j = 1:nSlices,
        [im1a_ph(:,:,:,j,k),centPos] = nnUnwrap(squeeze(im1a_r(:,:,j,:,k)),squeeze(im1a_i(:,:,j,:,k)),1,centPos,showFig);
    end
end

im_ph_P = im1a_ph(:,:,:,:,1)-im1a_ph(:,:,:,:,4);
im_ph_M = im1a_ph(:,:,:,:,2)-im1a_ph(:,:,:,:,4);
im_ph_S = im1a_ph(:,:,:,:,3)-im1a_ph(:,:,:,:,4);
mn_P = mean(im_ph_P,3);
mn_M = mean(im_ph_M,3);
mn_S = mean(im_ph_S,3);
% swap slice and phase directions and subtract off direction-averaged mean phase
for k = 1:nPhases,
    for sl = 1:nSlices,
        im_phase(:,:,sl,k,1) = im_ph_P(:,:,k,sl)-mn_P(:,:,1,sl);
        im_phase(:,:,sl,k,2) = im_ph_M(:,:,k,sl)-mn_M(:,:,1,sl);
        im_phase(:,:,sl,k,3) = im_ph_S(:,:,k,sl)-mn_S(:,:,1,sl);
    end
end

im_tot = sqrt(im_ph_P.^2 + im_ph_M.^2 + im_ph_S.^2);

for ph = 1:nPhases,
    for sl = 1:nSlices, 
        im_mag(:,:,sl,ph) = mean(im1a_r(:,:,sl,ph,:),5);
    end
end


