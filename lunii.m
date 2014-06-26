function [nii_file, fileLocation] = lunii_multi(openText,fileLocation)

%% Function to load nii file (untouch)

if isempty(fileLocation),
    [f p] = uigetfile('*.nii',openText);
    fileLocation = [p f];
end

nii_file = load_untouch_nii(fileLocation);
