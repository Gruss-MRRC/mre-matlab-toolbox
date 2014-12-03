function imgf = nthHarmonic3D(img,n)
% Isolate n-th harmonic of 4D matrix,
% where 4th axis is time
% Usage:
% imgf = nthHarmonic(img,n)
% img ..... 4D matrix, 3 spatial dimensions, 1 temporal dimension
% n ....... natural number, get Nth harmonic
% imgf .... 3D amplitude matrix
  timeaxis=4;
  timesteps=size(img,timeaxis);
  if timesteps>1,
    imgf=fft(img,timesteps,timeaxis);
    imgf=imgf(:,:,:,n);
  else
    imgf=img; % nothing to do, already 3D
  end
