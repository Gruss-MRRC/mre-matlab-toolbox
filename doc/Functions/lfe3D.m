function [elastogramVol,varargout] = lfe3D(px,py,pz,mask,f,fov)
% Apply LFE reconstruction to 3D MRE volumes
%
% Usage:
% [elastogramVol,ex,ey,ez] = lfe3D(px,py,pz,m,f,fov)
%
% px,py,pz ....... 3D or 4D matrices (X,Y,Z,[T]) containing unwrapped phase 
%                  data for each sensitization direction (x,y,z encoding)
% mask............ 3D logical matrix                   
% f .............. mechanical frequency (Hz)
% fov ............ field-of-view ([X Y]) (meters)
% 
% elastogramVol .. 3D matrix of averaged elasticity values (kPa)
% ez,ey,ez      .. [optional] 3D matrices of elasticity values in each direction
% 
% LFE is applied to 2D slices along the third dimension for each phase matrix. 
% The 3D volume is the average of the results for each phase encoding direction,
% masked by the provided magnitude image. 
% 
% Author:
% Alex Krolick <amk283@cornell.edu>
%
% See also: lfe

lfeOut = zeros(size(mask));
lfeOutX = lfeOut; lfeOutY = lfeOut; lfeOutZ = lfeOut;
nSlices = size(px,3);

for z=1:nSlices
  lfeOutX(:,:,z) = lfe(squeeze(px(:,:,z,:)),f,fov,16,3,2,2);
  lfeOutY(:,:,z) = lfe(squeeze(py(:,:,z,:)),f,fov,16,3,2,2);
  lfeOutZ(:,:,z) = lfe(squeeze(pz(:,:,z,:)),f,fov,16,3,2,2);
end

% Merge by averaging results for each direction
lfeOut = (lfeOutX+lfeOutY+lfeOutZ)/3;

% Mask
elastogramVol = lfeOut.*(mask>0).*(lfeOut<5);

% Parse optional outputs
nout = max(nargout,1)-1;
if nout == 3,
	varargout(1) = {lfeOutX};
	varargout(2) = {lfeOutY};
	varargout(3) = {lfeOutZ};
end

% TODO
% - Take only one input and check dimensions
% --> Allow MxNxSx1, MxNxSxP, MxNxSxPxV, MxNxSxPxV

