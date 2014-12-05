function amplitude = bulkmotion(P,scale,direction,showplots)
% Compute amplitude of bulk motion waveform from a phase map. Use the mean
% phases of image slices (3rd dimension) to infer the bulk motion waveform
% (4th dimension) of the entire volume. Find the peak-to-peak amplitude of
% this wave.
% 
% Usage:
% amplitude = bulkmotion(P,scale,direction,showplots)
%
% Inputs:
% P .......... Phase map (4 or 5D), where the dimensions are
%              (X,Y,Z,time,[direction]) and units are radians
% scale ...... Conversion factor between radians and output units
% direction .. If P is 5D, which direction to use
%              If P is 4D, must equal 1
% showplots .. Boolean; whether or not to open a figure with motion plots
%
% Outputs:
% amplitude .. Peak-to-peak wave amplitude (in radians) * scale


[lx,ly,lz,nt,nd] = size(P); % length in x,length in y,length in z, 
                            % number of time (phase) steps, 
                            % number of dimensions

% Get mean phases for each layer and timestep:
meanphases = zeros(lz,nt); 
for i=1:lz, 
  for j = 1:nt,
  % Region of interest (roi) defined as rectangle inset by a fraction of
  % side length on each side:
  inset = 0;
  if inset>0 && inset<0.5
	roi = P(ceil(lx*inset):ceil(lx-lx*inset),...
          ceil(ly*inset):ceil(ly-ly*inset),...
          i,j,direction);
  else 
  	roi=P(:,:,i,j,direction);
  end
  % Change to vector of unmasked elements only (in case NaN mask was used):
  roi = roi(roi~=0);
  roi = roi(~isnan(roi));
  meanphases(i,j) = mean(roi);
  end
end
% Unwrap (correct jumps of pi or more between time steps):
%meanphases = unwrap(meanphases);
% Take mean in time for each slice and subtract it to center each wave
% about 0. This corrects for phase offsets between layers.
temporalmeans = mean(meanphases');
for i = 1:lz, meanphases(i,:) = meanphases(i,:)-temporalmeans(i); end
% Unwrap the normalized mean-phase waveform for each slice:
%meanphases = unwrap(meanphases');
% Use the median wave as the bulk motion waveform for the entire brain:
bulkwave = median(meanphases);
% Show plots, if enabled
if showplots
  figure,hold on
  plot((meanphases'),':')
  plot(bulkwave,'r-','LineWidth',5)
end
% Calculate peak to peak motion amplitude
phaseamp = max(bulkwave)-min(bulkwave);
amplitude = phaseamp*scale;
