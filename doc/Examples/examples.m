%% MRE Examples
% Using the MRE Matlab Toolbox
%%

%% Open an MRI image
% MRI2MAT opens common image formats such as DICOM and NIFTI
%
% * <matlab:doc('mri2mat') Help for mri2mat.m>
% * <matlab:run('mri2mat') Run mri2mat.m>
%
%   >> [M,info] = mri2mat();
%
% <<filepicker.png>>
[M,info] = mri2mat();

%%
% The file you select in the file browser is opened as a matrix, which we
% can see has the same dimensions as those specified in the NIFTI header:
%
%   info.dime.dim(1)
%   length(size(M))
disp(info.dime.dim(1))
disp(length(size(M)))
%%
%   info.dime.dim(2:5)
%   size(M)
disp(info.dime.dim(2:5))
disp(size(M))


%% Visualization in 2D
% View a 2D slice of the image using the Image Toolbox function _imagesc_
%
%   im = flipud(squeeze(M(65,:,:,1))');
%   imagesc(im)
%   colormap gray, axis tight

  im = flipud(squeeze(M(65,:,:,1))');
  imagesc(im)
  colormap gray, axis tight
  
%%
% Note the use of _flipud_ and _transpose (')_ to orient the
% image and _squeeze_ to reduce the image to 2D


%% Visualization in 3D
% View a 3D volume
%
%   S = smooth3(M(:,:,:,1));
%   isosurface(S,0)
%
% <<isosurface.png>>

%% Section plane through a 3D volume
% Create a section plane through a 3D volume, using information from the 
% NIFTI header to scale the axes
%
%   colormap gray
%   D = M(round(1:info.dime.dim(2)/2),:,:,1); % split along sagittal plane
%   Ds = smooth3(D);
%   hiso = patch(isosurface(Ds,0),...
%     'FaceColor',[.5,.5,.8],...
%     'EdgeColor','none');
%   isonormals(Ds,hiso)
%   hcap = patch(isocaps(D,0),...
%     'FaceColor','interp',...
%     'EdgeColor','none');
%   view(35,30) 
%   axis tight 
%   pixdim = info.dime.pixdim(2:4); % get voxel dimensions
%   daspect(1./[pixdim(1),pixdim(2),pixdim(3)]) % scale proportional to voxel shape
%   lightangle(90,0);
%   set(gcf,'Renderer','opengl'); lighting phong
%   set(hcap,'AmbientStrength',1)
%   set(hiso,'SpecularColorReflectance',0,'SpecularExponent',50)
%
% <<brainsection.png>>

%% MREView
% Use _MREView_ to visualize volumes of up to 5 dimensions
%
% * <matlab:doc('MREView') Help for MREView.m>
% * <matlab:run('MREView') Run MREView.m>
%
%   >> MREView
%
% <<mreview.png>>
%
% The _unwrap_ button launches a phase unwrapping function that attempts to
% correct large jumps in the aparent phase due to aliasing.
%
% * <matlab:doc('fastUnwrap') Help for fastUnwrap>
% * <matlab:doc('QualityGuidedUnwrap2D_r1') Help for QualityGuidedUnwrap>
% 
% <<unwrapdialog.png>>
%
% <<mreviewunwrap.png>>


%% Volume segmentation
% _segmentVolumes_ can be used to isolate volumes (e.g., the lateral
% ventricles) inside an MRI file. It has a number of GUI tools to
% facilitate image processing, and can create 3D views on demand.
% 
% * <matlab:doc('segmentVolumes') Help for segmentVolumes.m>
% * <matlab:run('segmentVolumes') Run segmentVolumes.m>
% 
% *Image Processing tools available in _segmentVolumes_*
%
% * <matlab:doc('imerode') imerode>
% * <matlab:doc('imdilate') imdilate>
% * <matlab:doc('bwlabel') bwlabel>
%
% The current volume and file information can be saved to the spreadsheet
% specified in the textbox at the top of the window by clicking _Save_ in
% the toolbar. This also saves changes to the file in Matlab's |mat| format
% inside the image's parent folder. 
%
%   >> segmentVolumes
%
% <<segmentvolumes.png>>
%
% Using the isolate tool to select a region
%
% <<isolatetool.png>>
%
% The image after running the isolate tool
%
% <<segmentvolumesisolate.png>>
%
% 3D rendering of isolated ventricles
%
% <<isolatedventricles3d.png>>

%% Make a movie of a slice
% * <matlab:doc('mreslice2avi') Help for mreslice2avi.m>
%
% *Usage:*
%
%   mreslice2avi(phase_image,movie_filename_prefix,slice)
%   >> data = load('Datasets/imgXYZ.mat')
%   >> mreslice2avi(data.P,'imgXYZ_slice',15)
data = load('~/Datasets/MRE/104A_MRE_4D.mat')
mreslice2avi(data.P,'demo_slice',15);
