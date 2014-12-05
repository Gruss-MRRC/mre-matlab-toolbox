function compareTissue(grayIn,whiteIn,shearIn,grayOut,whiteOut)
% Compare shear stiffness of gray/white matter using segmented elastogram
%
% Usage:
% compareTissue(grayIn,whiteIn,shearIn,grayOut,whiteOut)
%
% Input arguments (filenames):
% gray/whiteIn ... anatomy tissue masks, 3D Nifti
% shearIn ........ elastogram, 3D volume in 2D delimited data file
%                  as output by MRE/Wave
% gray/whiteOut .. output filenames for gray/white elastograms
%
% All output is converted to 32-bit float datatype 
% (the 'single' class in Matlab)

graymatter=load_untouch_nii(grayIn);
whitematter=load_untouch_nii(whiteIn);
shearmap=load(shearIn);
shearmap=reshape(shearmap',size(graymatter.img));
graymatter.img=single(...
					shearmap.*double(graymatter.img>0.2));
whitematter.img=single(...
					 shearmap.*double(whitematter.img>0.2));
graymatter.hdr.dime.datatype=16; %float32 code in NIFTI spec
whitematter.hdr.dime.datatype=16;
save_untouch_nii(graymatter,grayOut);
save_untouch_nii(whitematter,whiteOut);
