function stepthru3d(matrix3d)
% Step through a 3d matrix slice by slice (xy planes, increment z)
% Holds colormap constant
% Press any key in the terminal to advance
% Usage:
% stepthru3d(matrix3d)
m=matrix3d;
n=size(m,3);
figure
cax = [min(m(:)),max(m(:))];
for ii=1:n
	imagesc(m(:,:,ii),cax);
	colorbar
	axis image;
	fprintf('\rslice %d/%d',ii,n)
	pause()
end
fprintf('\n')
