% GoldsteinUnwrap2D is a script to demonstrate the 2D Goldstein branch cut phase unwrapping algorithm.
%
%  Calls: PhaseResidues
%         BranchCuts
%         FloodFill
%
% References::
% 1. R. M. Goldstein, H. A. Zebken, and C. L. Werner, "Satellite radar interferometry:
%    Two-dimensional phase unwrapping," Radio Sci., vol. 23, no. 4, pp. 713-720, 1988.
% 2. D. C. Ghiglia and M. D. Pritt, Two-Dimensional Phase Unwrapping:
%    Theory, Algorithms and Software. New York: Wiley-Interscience, 1998.
%
% Inputs:  1. IM = Complex image
%  Optional inputs:
%          2. im_mask = Binary mask
%          3. max_box_radius= Maximum search box radius (pixels)
%          4. threshold_std = Number of noise standard deviations used for
%          thresholding the magnitude image
% Outputs: 1. Unwrapped phase image
%          2. Phase quality map
%
% This code can easily be extended for 3D phase unwrapping.
%
% Posted by Bruce Spottiswoode on 22 December 2008
%
% 2010/07/23  Modified by Carey Smith
%             1. Changed IM_mask, IM_mag, IM_phase to lowercase, in order to be
%                similar to the Quality Guided routines.
%             2. Moved the logic to chose the reference point from FloodFill.m
%                to here, in order to be similar to the Quality Guided routines.
%                The coordinates of the reference point, colref, rowref, are now
%                passed to FloodFill2.m
%             3. User can optionally have the code automatically chose the
%             largest magnatude point as the reference point.
%             4. Allowed the user to specify IM, im_mask, max_box_radius before
%                calling this routine.  (Note that threshold_std is not used.)
%             5. Modified the plots to suit my preferences.
% 2014/07/25  Modified by Alex Krolick:
%             1. Change to function instead of script
%             2. Instead of complex image, expect MAG and PHASE matrices
%             3. Turn off interactivity and print statements
%             

function im_unwrapped = GoldsteinUnwrap2D_r1(im_mag,im_phase)
IM = im_mag; % dummy variable just used for sizing matrices

%% Replace with your mask (if desired)
im_mask = ones(size(IM));
mag_max = max(im_mag(:));
indx1 = im_mag < 0.1*mag_max;  %Intensity = mag^2, so this = .01 threshold on the intensity
im_mask(indx1) = 0;   % Don't mask at this point; wait until residues are computed
if(~exist('im_mask','var'))
  im_mask = ones(size(IM));          %Mask (if applicable)
end
%figure; imagesc(im_mag.*im_mask),   colormap(gray), axis square, axis off, title('GS Initial masked magnitude'); colorbar;
%figure; imagesc(im_phase.*im_mask), colormap(gray), axis square, axis off, title('GS Initial masked phase'); colorbar;

%% Compute the residues
residue_charge = PhaseResidues_r1(im_phase, im_mask); % Calculate phase residues (Does not use mask)
%figure; imagesc(residue_charge), colormap(gray), axis square, axis off, title('GS Phase residues (charged)'); colorbar;

%% Compute the branch cuts
max_box_radius=floor(length(residue_charge)/2);  % Maximum search box radius (pixels)
%max_box_radius=4  % Maximum search box radius (pixels)
if(~exist('max_box_radius','var'))
  max_box_radius=4;  % Maximum search box radius (pixels)
end
% BranchCuts() ignores residues with mask == 0, so keep the entire mask == 1
branch_cuts = BranchCuts_r1(residue_charge, max_box_radius, im_mask); % Place branch cuts
%figure; imagesc(branch_cuts),    colormap(gray), axis square, axis off, title('GS Branch cuts'); colorbar;

im_mask(branch_cuts) = 0;  % Now need to mask off branch cut points, in order to avoid an error in FloodFill
im_mag1 = im_mag.*im_mask; % Mask off magnitude == 0 points, so that they are not chosen for the starting point

%% Manually (default) or automatically identify starting seed point 
try   % Chose starting point = max. intensity
  [r_dim, c_dim]=size(im_phase);
  im_mag1(1,:) = 0;                     %Set magnitude of border pixels to 0, so that they are not used for the reference
  im_mag1(r_dim,:) = 0;
  im_mag1(:,1) = 0;
  im_mag1(:,c_dim) = 0;
  [rowrefn,colrefn] = find(im_mag1 >= 0.99*mag_max);
  rowref = rowrefn(1);                 %Choose the 1st point for a reference (known good value)
  colref = colrefn(1);                  %Choose the 1st point for a reference (known good value)
  if (branch_cuts(rowref,colref)==1)
    error(['Selected point: (',int2str(rowref),',',int2str(colref),' corresponds to a branch cut.']);
  end
catch
  im_phase_quality = im_mag1;
  minp = im_phase_quality(2:end-1, 2:end-1); minp = min(minp(:));
  maxp = im_phase_quality(2:end-1, 2:end-1); maxp = max(maxp(:));
  %uiwait(msgbox('Select known true phase reference phase point. Black = high quality phase; white = low quality phase.','Phase reference point','modal'));
  uiwait(msgbox('Select known true phase reference phase point. White = high magnitude; Black = low magnitude.','Phase reference point','modal'));
  h = figure; imagesc(im_phase_quality,[minp maxp]), colormap(gray), colorbar, axis square, axis off; title('Phase quality map');
  [xpoint,ypoint] = ginput(h);        %Select starting point for the guided floodfill algorithm
  colref = round(xpoint);
  rowref = round(ypoint);
  close(h);     pause(.1)                         %Close the figure;
end

%% Unwrap
if(exist('rowref','var'))
  im_unwrapped = FloodFill_r1(im_phase, im_mag, branch_cuts, im_mask, colref, rowref); % Flood fill phase unwrapping
else
  im_unwrapped = FloodFill_r1(im_phase, im_mag, branch_cuts, im_mask); % Flood fill phase unwrapping
end
% Display results
%figure; imagesc(im_unwrapped), colormap(gray), colorbar, axis square, axis off, title('GS Unwrapped phase');
