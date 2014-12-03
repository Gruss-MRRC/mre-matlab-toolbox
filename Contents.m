% MR Elastography Toolbox
%
% A Matlab toolbox for working with magnetic resonance elastography (MRE) images
%
% Examples
%   examples - MRE toolbox usage examples
%
% Functions
%   bandpassfilter            - Creates bandpass Butterworth filter in two dimensions.
%   bulkmotion                - Compute amplitude of bulk motion waveform from a phase map.
%   compareTissue             - Compare shear stiffness of gray/white matter using segmented elastogram
%   fixPhase                  - 
%   getFFT                    - 
%   getMRE_EPI                - 
%   getMREimages              - Extract magnitude and phase images from DICOM files.
%   getMRESinkus              - Extract magnitude and phase images from DICOM files.
%   getUmap                   - 
%   lfe                       - Create elastogram from MRE phase data using local-frequency 
%   lfe3D                     - Apply LFE reconstruction to 3D MRE volumes
%   lowpassfilter             - lowpassfilter.m - Creates lowpass Butterworth filter in two dimensions.
%   lunii                     - Function to load nii file (without applying affine xform)
%   mreslice2avi              - Save MRE motion sequence to avi movie.
%   mri2mat                   - Convert NIFTI and DICOM images to matrix format.
%   nthHarmonic               - Isolate n-th harmonic of phase map
%   nthHarmonic3D             - Isolate n-th harmonic of 4D matrix,
%   probePhaseData            - b0valsPos = [1 8 18 29 43 54 65 76 87 98];
%   projectMagnitudeOntoWalls - Make the "walls" of the 3D plot area show midvolume axial, coronal,
%   scannerphase2rad          - Convert units of phase image from MRI scanner (range 0...2^n) to radians
%   splitDimension            - Split matrix along dimension
%   stepthru3d                - Step through a 3d matrix slice by slice (xy planes, increment z)
%   unwrapAndNormalize        - Unwrap with ImageJ Phase Tools (interactive)
%   visualizeFlow             - Visualize out-of-plane fluid flow as a moving 3D surface
%   visualizeMotion           - Visualize 3D motion of a matrix slice
%   visualizeMotionTriPlanar  - Visualize motion of a matrix slice in 3D
%   Unwrap/*                  - Phase unwrapping utilities (some 3rd party under BSD license)
%   Plot/*                    - Plotting utilities (some 3rd party under BSD license)
%
% GUIs
%   MREBulkMotionCalculator - Calculate bulk motion amplitude of a cerebral MRE sequence
%   MREView                 - Preview and phase unwrap MRE sequences.
%   segmentVolumes          - Isolate and calculate volume of structures in MRI images.
%
% Scripts
%   MRE_sequence_simulation - Simulate an MRE acquisition sequence
%   prepareImages.sh        - Shell pipeline to preprocess MRE images and segment using additional T1 images
%   prepareImagesT2.sh      - Shell pipeline to preprocess MRE images and segment using T2 magnitude volumes
