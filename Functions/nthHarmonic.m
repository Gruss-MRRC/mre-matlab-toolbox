function imgf = nthHarmonic(img,n)
% Isolate n-th harmonic of phase map
% Usage:
% imgf = nthHarmonic(img,n)
% img ..... 3D matrix, 2 spatial dimensions, 1 temporal dimension
% n ....... natural number, get Nth harmonic
% imgf .... 2D amplitude matrix

  timeaxis=3;
  timesteps=size(img,timeaxis);
  if timesteps>1,
    imgf=fft(img,timesteps,timeaxis);
    imgf=imgf(:,:,n);
  else
    imgf=img; % nothing to do, already 2D
  end
  
  % TODO typecheck inputs
