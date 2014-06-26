function [im_mag, im_phase] = getMRESinkus(showFig)

[f p] = uigetfile('*.dicom');
filename = [p '../RAW/' strtok(f,'.') '.nii'];
if exist(filename) > 0,
    im1 = lunii('Select nifti image',filename);
else
    filename = [p '../RAW/' strtok(f,'.') '.nii.gz'];
    im1 = lunii('Select nifti image',filename);
end
im1 = im1.img;

nSlices = size(im1,3);
[status,result] = system(['dcmdump ' p f ' | grep NumberOfTemporalPositions | awk ''{print $3}'' | head -n 1 | sed ''s/\[//g'' | sed ''s/]//g''']);
nDirs = str2num(result);

nPhases = size(im1,4) / 2 / nDirs;

brainFilename = [p '../RAW/' strtok(f,'.') 'Brain.nii.gz'];
if exist(brainFilename) > 0,
    imBrain = lunii('Select NIFTI BET image',brainFilename);
else
    imBrain = lunii('Select NIFTI BET image','');
end
imBrain = imBrain.img;

k1=1;
for dir = 1:nDirs,
    for ph = 1:nPhases,
        im1a_r(:,:,:,ph,dir) = im1(:,:,:,k1) .* int16(imBrain>0);
        im1a_i(:,:,:,ph,dir) = im1(:,:,:,size(im1,4) / 2 + k1) .* int16(imBrain>0);
        k1 = k1 + 1;
    end
end

centPos = [];
for k = 1:nDirs, 
    for j = 1:nSlices,
        [im1a_ph(:,:,:,j,k),centPos] = nnUnwrap(squeeze(im1a_r(:,:,j,:,k)),squeeze(im1a_i(:,:,j,:,k)),1,centPos,showFig);
    end
end

% permute, so that slices are first
im1a_ph = permute(im1a_ph,[1 2 4 3 5]);

im_ph_P = im1a_ph(:,:,:,:,1)-im1a_ph(:,:,:,:,4);
im_ph_M = im1a_ph(:,:,:,:,2)-im1a_ph(:,:,:,:,4);
im_ph_S = im1a_ph(:,:,:,:,3)-im1a_ph(:,:,:,:,4);
mn_P = mean(im_ph_P,4);
mn_M = mean(im_ph_M,4);
mn_S = mean(im_ph_S,4);

for h = 1:8,
    mn_P(h) = mean(nonzeros(im_ph_P(:,:,4,h)));
    im_ph_P(:,:,:,h) = im_ph_P(:,:,:,h) - mn_P(h);
    mn_M(h) = mean(nonzeros(im_ph_M(:,:,4,h)));
    im_ph_M(:,:,:,h) = im_ph_M(:,:,:,h) - mn_M(h);
    mn_S(h) = mean(nonzeros(im_ph_S(:,:,4,h)));
    im_ph_S(:,:,:,h) = im_ph_S(:,:,:,h) - mn_S(h);
end

% subtract off direction-averaged mean phase
for ph = 1:nPhases,
    im_phase(:,:,:,ph,1) = im_ph_P(:,:,:,ph)-mn_P;
    im_phase(:,:,:,ph,2) = im_ph_M(:,:,:,ph)-mn_M;
    im_phase(:,:,:,ph,3) = im_ph_S(:,:,:,ph)-mn_S;
end

im_tot = sqrt(im_ph_P.^2 + im_ph_M.^2 + im_ph_S.^2);

for ph = 1:nPhases,
    for sl = 1:nSlices, 
        im_mag(:,:,sl,ph) = mean(im1a_r(:,:,sl,ph,:),5);
    end
end


