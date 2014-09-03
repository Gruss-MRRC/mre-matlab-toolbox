function mreslice2avi(phase_image,movie_filename_prefix,slice)

% Save MRE motion sequence to avi movie.
%
% Plot X,Y,Z components of motion side-by-side, capture each frame, and
% save to a movie file. Assumes a 5D phase image, where the dimensions are
% (x,y,slice,phase,direction). Phase is equivalent to timestep. Each frame
% in the movie is a step along the phase axis. 
%
% Inputs
% ------
% phase_image: 5D phase image
% movie_filename_prefix: String to prepend to the slice number in .avi file
% slice: Integer index of the slice (z-coordinate) to be viewed.
%        If you don't know, try ceil(size(phase_image,3)/2)
%
% Notes on creating 'phase_image'
% -------------------------------
% The MREView GUI can be used to preview DICOM and NIFTI files and launch
% phase unwrapping sequences for aliased images. Use the 'Save' button to
% store the processed image in Matlab's .mat format, and the 'load' command
% to pull in the variables for use in MRESLICE2AVI.
% 
% Authors
% -------
% Alex Krolick <amk283@cornell.edu>
%
% See also MREView

% Configuration
p = phase_image; % your phase image here
moviefilename = [movie_filename_prefix num2str(slice) '.avi'];
nphases = size(P,4); % number of phases (timesteps)

h = figure;
colormap gray
for i = 1:nphases
  subplot(1,3,1)
  imagesc(p(:,:,slice,i,1))
  axis square
  axis off
  title('X')
  subplot(1,3,2)
  imagesc(p(:,:,slice,i,2))
  axis  square
  axis off
  title('Y')
  subplot(1,3,3)
  imagesc(p(:,:,slice,i,3))
  axis square
  axis off
  title('Z')
  F(i) = getframe(h);
end
movie2avi(F, moviefilename);
