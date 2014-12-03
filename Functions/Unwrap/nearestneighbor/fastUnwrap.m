function P = fastUnwrap(M,P)
% Unwrap the phase of a multidimensional image. 
%
% Matlab's `unwrap` function operates columnwise and doesn't work well
% for data >1D. This works on 2D to 5D images and, in particular, is
% designed to unwrap magnetic resonance elastography (MRE) images stored in
% DICOM or NIFTI format.
%
% Usage
% ----
% P = fastUnwrap(M,P)
%
% Inputs
% ------
% M: magnitude image
% P: phase image (-pi,pi)
%
% If either M or P are empty ([]), create using `mri2mat`. M and P can have
% up to 5 dimensions. If P is being unwrapped and has the dimensions
% P(x,y,slice,phase,direction), the algorithm iterates over slice, phase,
% and direction while unwrapping the xy plane.
%
% Outputs
% -------
% P: Unwrapped phase image
%
% Algorithm
% ---------
% Nearest-neighbor unwrapping using a block of configurable size that
% follows a spiral path. Inside the block, +/- 2pi is added to elements
% that deviate by more than pi from the central value. The path generation
% and unwrapping algorithm are subfunctions that can easily be substituted
% for other methods.
%
% Limitations
% -----------
% This is a very naive algorithm that breaks down on noisy or heavily wrapped
% data. More robust, but slower algorithms are provided for working with these
% datasets in Matlab. If speed is a priority, a Mex call to a C library or use
% of a dedicated image processor such as ImageJ is recommended.
%
% Author
% ------
% Alex Krolick <amk283@cornell.edu>
%
% See also mri2mat, MRE_Preview, unwrap

  %--Config--%
  verbose=false; % Toggle debugging flags
  
  % If M or P are empty, open a new file using mri2mat
  if (isempty(P) || isempty(M)), [M,P] = mri2mat(); end
  
  [nX,nY,nSlices,nPhases,nDirs] = size(P);
  %--Variables to control path--%
  blockSize = 3;
  step      = floor(blockSize/2);
  center    = [round(nX/2),round(nY/2)];
  [X,Y]     = pathGen(center,1,nX-blockSize);
  
  %--Mask parameters--%
  % Cut off everything below 12.5% of peak magnitude
  method='minimum';
  threshold=.125*(max(M(:))-min(M(:)));
  maskvalue=NaN;
  
  %--Unwrap along path--%
  for dir = 1:nDirs, for slice = 1:nSlices, for phase = 1:nPhases,
%     Ps = mask(... % mask current slice with magnitude filter
%         P(:,:,slice,phase,dir),...
%         M(:,:,slice,phase,dir),...
%         method,threshold,...
%         maskvalue);
%     Ps(Ps==-pi)=NaN;
    Ps=P(:,:,slice,phase,dir);
    x=X(1); y=Y(1);
    i=1;
    while i<length(X)
      xb = (x-step):(x+step);
      yb = (y-step):(y+step);
      Ps(xb,yb) = unwrapBlock(Ps(xb,yb),Ps(x,y));
      i=i+1; x=X(i); y=Y(i);
    end
    if(verbose), debug(Ps), end
    tmp=Ps(~isnan(Ps));
    mn=mean(nonzeros(tmp(:)));
    P(:,:,slice,phase,dir) = Ps-mn;
  end, end, end % end for-loops

  
function A = mask(A,B,mode,threshold,value)
% Apply a mask to `A` if the corresponding value in `B` is below the
% threshold parameter. `Value` is the new value that maskeded-out elements
% will take. A and B must be the same size.
% Modes (element in A gets `value` if...):
% factor:       element in B less than `threshold` times B's peak value
% minumum:      element in B less than `threshold`
% rangeExclude: element in B outside of `threshold` range (threshold is a
%               2-element vector)
% rangeInclude: element in B inside of `threshold` range (threshold is a
%               2-element vector)
% stdDev:       element in B more than `threshold` standard deviations
%               below the mean in B
% binary:       element in B equals 0 (B is a binary mask)

  switch mode
    case 'factor'
      m = max(B(:));
      msk = B<threshold*m;
    case 'minimum'
      msk = B<threshold;
    case 'rangeExclude'
      msk = B>threshold(1) & B<threshold(2);
    case 'rangeInclude'
      msk = B<threshold(1) | B>threshold(2);
    case 'stdDev'
      msk = B<(mean(nonzeros(B))-threshold*std(nonzeros(B)));
    case 'binary'
        msk = ~B;
  end
  msk = ~imfill(~msk,conndef(length(size(A)),'maximal'),'holes');
  % NOTE Mask elements are 1 if masked, but imfill thinks 0s are background
  % NOTE conndef updates the connectivity required to define a hole based
  % on the dimension of the matrix (2D is 9-1=8, 3D is 27-1=26, etc.)
 
  A(msk) = value;
  
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
      while abs(diff)>pi && tries<2
      diff2 = center-block(i)-sign(diff)*2*pi;
      if abs(diff2)<abs(diff)
        block(i) = block(i)+sign(diff)*2*pi;
      end
      diff = center-block(i);
      tries=tries+1;
      end
    end
  end
  
 
function block = moduloBlock(block,center)
% Perform phase unwrapping step on the provided subset (`block`) of the
% whole image
block = mod(block,2*pi);
  
  
function debug(Ps)
% Display verbose output, calculate gradients, show progress of unwrapping
% in graphics window
  [Py,Px] = gradient(Ps);
  P_ = sqrt(Px.^2+Py.^2);
  errors = numel(nonzeros(abs(P_)>pi));
  disp(errors)
  show(Ps)

