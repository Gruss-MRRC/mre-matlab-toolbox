% GuidedFloodFill.m unwraps a single 2D image using the quaility guided technique
% function IM_unwrapped = GuidedFloodFill(IM_phase, IM_mag, IM_unwrapped, unwrapped_binary, derivative_variance, adjoin, IM_mask)
% It can also be used to unwrap phases in cases when there are no branch cuts.
%
% Input: IM_phase, IM_unwrapped (seed points / pixels already unwrapped),
% unwrapped_binary the derivative variance, an adjoining matrix and a mask.
% Inputs:
%  IM_phase     = 2D array of wrapped phases (rads)
%  IM_mag       = 2D array of magnitudes
%  IM_unwrapped = a 2D array for the unwrapped phases (rads)
%                 This is initialized for the reference point set.
%  unwrapped_binary = a 2d array to identifying the points that have been unwrapped.
%                 This is initialized for the reference point set.
%  derivative_variance = a 2D array of the IM_phase's derivative variances
%  adjoin       = a 2D identifying the valid pixels adjoining the unwrapped pixels.
%                 This is initialized for the reference point set.
%  IM_mask      = a 2D array identifying the valid pixels
% Outputs:
%  IM_unwrapped = the 2D array of unwrapped phases (rads)
%
% This code can easily be extended for 3D phase unwrapping.
% Technique adapted from:
% D. C. Ghiglia and M. D. Pritt, Two-Dimensional Phase Unwrapping:
% Theory, Algorithms and Software. New York: Wiley-Interscience, 1998.

% Created by B.S. Spottiswoode on 11/11/2004
%
% 2010/07/23  Carey Smith
%             1. Implemented Itoh's in-lin method to remove 2*pi jumps
%                (rather than calling unwrap).
%             2. When there are 2 or more valid neighbors, use the one with the 
%                largest magnitude for determining the unwrapped phase.
%                (Previous code re-computed the unwrapped phase for each valid 
%                neighbor, over-writing the previous computation each time.)
%                The new logic eliminated the need for isolated_adjoining_pixel_flag.
%             3. Allow the edge pixels to be unwrapped.
%             4. Implemented Content Advisor's recommended usage of & vs. &&
%                and | vs. ||
% 2010/07/26  Carey Smith: Made a speed improvement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function IM_unwrapped = GuidedFloodFill2(IM_phase, IM_mag, IM_unwrapped, unwrapped_binary, derivative_variance, adjoin, IM_mask)

[r_dim, c_dim] = size(IM_phase);

% Include edge pixels
while sum(sum(adjoin(:))) ~= 0  %Loop until there are no more adjoining pixels
  adjoining_derivative_variance = derivative_variance.*adjoin + 101.*~adjoin; %Derivative variance values of the adjoining pixels (pad the zero adjoining values with 100)
  min_deriv_var = min(adjoining_derivative_variance(:)); % the minimum derivative variance
  if(min_deriv_var >= 101)  % 101 is an indicator of already assigned variance
    break;  % finished
  end
  [r_adjoin, c_adjoin] = find(adjoining_derivative_variance==min_deriv_var); %Obtain coordinates of the adjoining unwrapped phase pixel w/ the minimum derivative variance
  for(ii = 1:length(r_adjoin))
    r_active = r_adjoin(ii);
    c_active = c_adjoin(ii);
    if(adjoin(r_active,c_active)==1)  % find a valid point
      break;
    end
  end
  if(adjoin(r_active,c_active)==0)  % No valid point was found, 
    break;
  end
  phasev   = nan(1,4);     % Initialize.  Will overwrite for valid pixels
  IM_magv  = nan(1,4);     % Initialize.  Will overwrite for valid pixels
  %First search below for an adjoining unwrapped phase pixel
  if(r_active+1<=r_dim)  % Is this a valid index?
    if unwrapped_binary(r_active+1, c_active)==1
      phase_ref = IM_unwrapped(r_active+1, c_active);       % Obtain the reference unwrapped phase
      % Itoh's Method (suggested by 'Eric' on MATLAB Central to use for a length 2 vector):
      D = IM_phase(r_active, c_active)-phase_ref;
      deltap = atan2(sin(D),cos(D));   % Make it modulo +/-pi
      phasev(1) = phase_ref + deltap;  % This is the unwrapped phase
      IM_magv(1)= IM_mag(r_active+1, c_active);
    else % unwrapped_binary(r_active+1, c_active)==0
      if(IM_mask(r_active+1, c_active)==1)
        adjoin(r_active+1, c_active) = 1;  % Put the elgible, still-wrapped neighbors of this pixels in the adjoin set
      end
    end
  end
  %Then search above
  if(r_active-1>=1)  % Is this a valid index?
    if unwrapped_binary(r_active-1, c_active)==1
      phase_ref = IM_unwrapped(r_active-1, c_active);                                   %Obtain the reference unwrapped phase
      D = IM_phase(r_active, c_active)-phase_ref;
      deltap = atan2(sin(D),cos(D));   % Make it modulo +/-pi
      phasev(2) = phase_ref + deltap;  % This is the unwrapped phase
      IM_magv(2)= IM_mag(r_active-1, c_active);
    else % unwrapped_binary(r_active-1, c_active)==0
      if(IM_mask(r_active-1, c_active)==1)
        adjoin(r_active-1, c_active) = 1;  % Put the elgible, still-wrapped neighbors of this pixels in the adjoin set
      end
    end
  end
  %Then search on the right
  if(c_active+1<=c_dim)  % Is this a valid index?
    if unwrapped_binary(r_active, c_active+1)==1
      phase_ref = IM_unwrapped(r_active, c_active+1);                                   %Obtain the reference unwrapped phase
      D = IM_phase(r_active, c_active)-phase_ref;
      deltap = atan2(sin(D),cos(D));   % Make it modulo +/-pi
      phasev(3) = phase_ref + deltap;  % This is the unwrapped phase
      IM_magv(3)= IM_mag(r_active, c_active+1);
    else % unwrapped_binary(r_active, c_active+1)==0
      if(IM_mask(r_active, c_active+1)==1)
        adjoin(r_active, c_active+1) = 1;  % Put the elgible, still-wrapped neighbors of this pixels in the adjoin set
      end
    end
  end
  %Finally search on the left
  if(c_active-1>=1)  % Is this a valid index?
    if unwrapped_binary(r_active, c_active-1)==1
      phase_ref = IM_unwrapped(r_active, c_active-1);                                   %Obtain the reference unwrapped phase
      D = IM_phase(r_active, c_active)-phase_ref;
      deltap = atan2(sin(D),cos(D));   % Make it modulo +/-pi
      phasev(4) = phase_ref + deltap;  % This is the unwrapped phase
      IM_magv(4)= IM_mag(r_active, c_active-1);
    else % unwrapped_binary(r_active, c_active-1)==0
      if(IM_mask(r_active, c_active-1)==1)
        adjoin(r_active, c_active-1) = 1;  % Put the elgible, still-wrapped neighbors of this pixels in the adjoin set
      end
    end
  end
  idx_del = ~isnan(phasev);
  if(any(idx_del)) % Any valid adjoining pixels?
       % Use the strongest neighbor
    IM_max  = max(IM_magv(idx_del));
    idx_max = find((IM_magv >= 0.99*IM_max) & (idx_del==1));
    IM_unwrapped(r_active, c_active) = phasev(idx_max(1));  % Use the first, if there is a tie
    unwrapped_binary(r_active, c_active) = 1;      %Mark the pixel as unwrapped
    adjoin(r_active, c_active) = 0;              %Remove it from the list of adjoining pixels
  else  % no valid adjoining pixels found
    adjoin(r_active,c_active) = 0;  %Remove the current active pixel from the adjoin list
    continue;
  end
  %end
end % while sum(sum(adjoin(2:r_dim-1,2:c_dim-1))) ~= 0  %Loop until there are no more adjoining pixels
%disp(['All of the valid interior pixels have been calculated']);
