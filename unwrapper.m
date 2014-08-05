function P = unwrapper(M,P)
% Unwrap the phase of a multidimensional image. Matlab's `unwrap` function
% operates columnwise and doesn't work well for data >1D. This works on 2D
% to 5D images and, in particular, is designed to unwrap magnetic resonance
% elastography (MRE) images stored in DICOM or NIFTI format.
%
% Inputs:
% M:  magnitude image
% P:  phase image (-pi,pi)
%
% If either M or P are empty ([]), create using mri2mat. M and P can have
% up to 5 dimensions. If P is being unwrapped and has the dimensions
% P(x,y,slice,phase,direction), the algorithm iterates over slice, phase,
% and direction while unwrapping the xy plane.
%
% Outputs:
% P:  Unwrapped phase image
%
% Algorithm: 
% Nearest-neighbor unwrapping using a block of configurable size that
% follows a spiral path. Inside the block, +/- 2pi is added to elements
% that deviate by more than pi from the central value. The path generation
% and unwrapping algorithm are subfunctions that can easily be substituted
% for other methods.
%
% See also mri2mat, MRE_Preview, unwrap
 
  %--Config--%
  verbose=0; % Toggle debugging flags
  
  % If M or P are empty, open a new file using mri2mat
  if (isempty(P) || isempty(M)), [M,P] = mri2mat(); end
  
  P = mask(P,M,'minimum',180,NaN);
  [nX,nY,nSlices,nPhases,nDirs] = size(P);
  
  %--Variables to control path--%
  blockSize = 3;
  step      = floor(blockSize/2);
  center    = [round(nX/2),round(nY/2)];
  [X,Y]     = pathGen(center,1,nX-blockSize);
  
  %--Unwrap along path--%
  for dir = 1:nDirs, for slice = 1:nSlices, for phase = 1:nPhases,
    Ps = P(:,:,slice,phase,dir);
    x=X(1); y=Y(1);
    i=1;
    while i<length(X)
      xb = (x-step):(x+step);
      yb = (y-step):(y+step);
      Ps(xb,yb) = unwrapBlock(Ps(xb,yb),Ps(x,y));
      i=i+1; x=X(i); y=Y(i);
    end
    if(verbose), debug(Ps), end
    P(:,:,phase,slice,dir) = Ps;
  end, end, end % end for-loops

  
function A = mask(A,B,mode,threshold,value)
% Apply a mask to `A` if the corresponding value in `B` is below the
% threshold value in `B`. The threshold can be set to either a range, a
% percentage of the peak value, or a hard limit by setting `mode` to 
% 'range', 'factor', or 'minimum', respectively. `Value` is the value that 
% masked elements will take.
% A and B must be the same size
  switch mode
    case 'factor'
      m = max(B(:));
      A(B<threshold*m) = value;
    case 'minimum'
      A(B<threshold) = value;
    case 'rangeExclude'
      A(B>threshold(1) & B<threshold(2)) = value;
    case 'rangeInclude'
      A(B<threshold(1) | B>threshold(2)) = value;
    case 'stdDev'
      A(B<(mean(B)-threshold*std(B))) = value;
  end

  
function show(M)
% Display an image with display settings that make sense for MRI
  figure
  imagesc(M)
  axis equal
  colormap gray
  
  
function [x,y] = pathGen(center,scale,sideLength)
% Calculate x,y coordinates of a square spiral path with a space of size
% `scale` between loops
  n = 0;
  x = [center(1)]; y = [center(2)];
  target = 1;
  direction = 1;
  while y(end)<sideLength,
    while n<target,
      x= [x x(end)+direction*scale];
      y= [y y(end)];
      n= n+1;
    end
    n=0;
    while n<target,
      y= [y y(end)+direction*scale];
      x= [x x(end)];
      n= n+1;
    end
    n=0;
    target=target+1;
    direction = -direction;
  end


function block = unwrapBlock(block,center)
% Perform phase unwrapping step on the provided subset (`block`) of the
% whole image
  for i = 1:numel(block);
    if ~isnan(block(i))
      tries=1;
      diff = center-block(i);
      while diff>pi && tries<4
      diff2 = center-block(i)-sign(diff)*2*pi;
      if abs(diff)>abs(diff2)
        block(i) = block(i)+2*pi;
      end
      diff = center-block(i);
      tries=tries+1;
      end
    end
  end
  
  
function debug(Ps)
% Display verbose output, calculate gradients, show progress of unwrapping
% in graphics window
  [Py,Px] = gradient(Ps);
  P_ = sqrt(Px.^2+Py.^2);
  errors = numel(nonzeros(abs(P_)>pi));
  disp(errors)
  show(Ps)

