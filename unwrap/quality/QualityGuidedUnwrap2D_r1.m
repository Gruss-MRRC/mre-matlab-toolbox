%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QualityGuidedUnwrap2D implements 2D quality guided path following phase
% unwrapping algorithm.
%
% Inputs: 1. Complex image in .mat double format
%         2. Binary mask (optional)          
% Outputs: 1. Unwrapped phase image
%          2. Phase quality map
%
% This code can easily be extended for 3D phase unwrapping.
% Technique adapted from:
% D. C. Ghiglia and M. D. Pritt, Two-Dimensional Phase Unwrapping:
% Theory, Algorithms and Software. New York: Wiley-Interscience, 1998.
%
% Posted by Bruce Spottiswoode on 22 December 2008
% 2010/07/23  Modified by Carey Smith
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Replace with your images

% if(~exist('IM'))
%   load 'IM.mat'                               %Load complex image
% end
% im_mag   = abs(IM);                  %Magnitude image
% im_phase = angle(IM);                %Phase image

function im_unwrapped = QualityGuidedUnwrap2D_r1(im_mag,im_phase)
IM = ones(size(im_mag)); % Not really the complex image; just used for size()

%% Replace with your mask (if required)
mag_max = max(im_mag(:));
indx1 = im_mag < 0.1*mag_max;  %Intensity = mag^2, so this = .01 threshold on the intensity
%mag_std = std(im_mag(:)) % Use if you want standard-deviation-based mask
%mag_avg = mean(im_mag(:))
%indx1 = im_mag < (mag_avg-2*mag_std);  %2 std_devs below mean
im_mask = ones(size(IM));
im_mask(indx1) = 0;                  %Mask
if(~exist('im_mask','var'))
  im_mask = (ones(size(IM)));          %Mask (if applicable)
end
%figure; imagesc(im_mag.*im_mask),   colormap(gray), axis square, axis off, title('Initial masked magnitude'); colorbar;
%figure; imagesc(im_phase.*im_mask), colormap(gray), axis square, axis off, title('Initial masked phase'); colorbar;

im_unwrapped = nan(size(IM));        %Initialze the output unwrapped version of the phase
adjoin = zeros(size(IM));            %Zero starting matrix for adjoin matrix
unwrapped_binary = zeros(size(IM));  %Binary image to mark unwrapped pixels

%% Calculate phase quality map
im_phase_quality = PhaseDerivativeVariance_r1(im_phase);

%% Automatically (default) or manually identify starting seed point on a phase quality map 
minp = im_phase_quality(2:end-1, 2:end-1); minp = min(minp(:));
maxp = im_phase_quality(2:end-1, 2:end-1); maxp = max(maxp(:));
if(0)    % Chose starting point interactively
  figure; imagesc(im_phase_quality,[minp maxp]), colormap(gray), colorbar, axis square, axis off; title('Phase quality map'); 
  uiwait(msgbox('Select known true phase reference phase point. Black = high quality phase; white = low quality phase.','Phase reference point','modal'));
  [xpoint,ypoint] = ginput(1);                %Select starting point for the guided floodfill algorithm
  colref = round(xpoint);
  rowref = round(ypoint);
  close;  % close the figure;
else   % Chose starting point = max. intensity, but avoid an edge pixel
  [rowrefn,colrefn] = find(im_mag(2:end-1, 2:end-1) >= 0.99*mag_max);
  rowref = rowrefn(1)+1; % choose the 1st point for a reference (known good value)
  colref = colrefn(1)+1; % choose the 1st point for a reference (known good value)
end

%% Unwrap
im_unwrapped(rowref,colref) = im_phase(rowref,colref);                          %Save the unwrapped values
unwrapped_binary(rowref,colref,1) = 1;
if im_mask(rowref-1, colref, 1)==1;  adjoin(rowref-1, colref, 1) = 1; end       %Mark the pixels adjoining the selected point
if im_mask(rowref+1, colref, 1)==1;  adjoin(rowref+1, colref, 1) = 1; end
if im_mask(rowref, colref-1, 1)==1;  adjoin(rowref, colref-1, 1) = 1; end
if im_mask(rowref, colref+1, 1)==1;  adjoin(rowref, colref+1, 1) = 1; end
im_unwrapped = GuidedFloodFill_r1(im_phase, im_mag, im_unwrapped, unwrapped_binary, im_phase_quality, adjoin, im_mask);    %Unwrap

%% Plot images
%figure; imagesc(im_mag),       colormap(gray), colorbar, axis square, axis off; title('QG Magnitude image'); 
%figure; imagesc(im_phase),     colormap(gray), colorbar, axis square, axis off; title('QG Wrapped phase'); 
%figure; imagesc(im_phase_quality,[minp maxp]), colormap(gray), axis square, axis off, title('QG Phase quality map'); colorbar;
%figure; imagesc(im_unwrapped), colormap(gray), colorbar, axis square, axis off; title('QG Unwrapped phase'); 
