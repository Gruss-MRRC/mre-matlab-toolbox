function magnitudeWalls(axHandle,magnitudeMat,transparentBG,sliceN)
% Make the "walls" of the 3D plot area show midvolume axial, coronal,
% and sagittal projections of an image
%
% Usage:
% magnitudeWalls(axHandle,magnitudeMat,transparentBG,sliceN)
% 
% axHandle .....  Handle of the axes to draw the walls 
%                 (can get this for the current axis with 'gca' command)
% magnitudeMat .. 3D matrix
% transparentBG . Boolean. Toggles whether to make the background 0s transparent
%                 Note: this makes the program a lot slower
% sliceN ........ Number. If not empty ([]), show this slice number instead 
%                 of the midvolume axial (XY) plane
%
% Author:
% Alex Krolick <amk283@cornell.edu>

h = axHandle;
M = magnitudeMat;
mid = ceil(size(M)/2);

if ~isempty(sliceN) && isscalar(sliceN)
	mid(3)=sliceN;
end

mm = slice(M,[mid(1)],[],[]);
nn = slice(M,[],[mid(2)],[]);
pp = slice(M,[],[],[mid(3)]);
set(mm,'XData',get(mm,'XData')-mid(1)+1);
set(mm,'EdgeColor','none');
set(nn,'YData',get(nn,'YData')-mid(2)+1);
set(nn,'EdgeColor','none');
set(pp,'ZData',get(pp,'ZData')-mid(3)+1);
set(pp,'EdgeColor','none');

if transparentBG
	set(nn,'AlphaData',double(get(nn,'CData')>100));
	set(nn,'FaceAlpha','flat')
	set(mm,'AlphaData',double(get(mm,'CData')>100));
	set(mm,'FaceAlpha','flat')
	set(pp,'AlphaData',double(get(pp,'CData')>100));
	set(pp,'FaceAlpha','flat')
end
