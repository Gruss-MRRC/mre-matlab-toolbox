function amplitude = bulkmotion(P,scale,direction,showplots)
% Compute amplitude of bulk motion waveform from a phase map. Use the mean
% phases of image slices (3rd dimension) to infer the bulk motion waveform
% (4th dimension) of the entire volume. Find the peak-to-peak amplitude of
% this wave.
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


[lx,ly,lz,nt,nd] = size(P); % x,y,slice,time (phase), dimensions

meanphases = zeros(lz,nt);
for i=1:lz, 
  for j = 1:nt,
  % Region of interest (roi) defined as rectangle inset by a fraction of
  % side length on each side:
  inset = 1/3; 
  roi = P(round(lx*inset):round(lx-lx*inset),...
          round(ly*inset):round(ly-ly*inset),...
          i,j,direction);
  % Change to vector of unmasked elements only (in case NaN mask was used):
  roi = roi(~isnan(roi));
  meanphases(i,j) = mean(roi);
  end
end
meanphases = unwrap(meanphases);
% Take mean in time for each slice and subtract it to center each wave
% about 0. This corrects for phase offsets between layers that may be due
% to a prior unwrapping step or a property of the acquistion:
temporalmeans = mean(meanphases');
for i = 1:lz, meanphases(i,:) = meanphases(i,:)-temporalmeans(i); end
% Unwrap the normalized mean-phase waveform for each slice:
meanphases = unwrap(meanphases');
% Use the median wave as the bulk motion waveform for the entire brain:
bulkwave = median(meanphases');
% Show plots, if enabled
if showplots
  figure,hold on
  plot((meanphases),':')
  plot(bulkwave,'r-','LineWidth',5)
end
% Calculate peak to peak motion amplitude
phaseamp = max(bulkwave)-min(bulkwave);
amplitude = phaseamp*scale;
