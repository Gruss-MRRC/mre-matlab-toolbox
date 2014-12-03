function unwrapAndNormalize(magfile,rawmagfile,phasefile,outputfile)
% 1. Unwrap with ImageJ Phase Tools (interactive)
% 2. Normalize by subtracting mean phase from each slice
% Inputs:
% magfile ..... NIFTI magnitude, brain extracted 3D volume.
%               Filename (string).
% rawmagfile .. NIFTI magnitude, 4D volume.
%               Filename (string).
% phasefile ... NIFTI phase, 4D volume sequence w/5th dim 
%               [sensitization direction] flattened and 
%               4 motion encoding gradients.
%               Units should be milliradians.
%               Filename (string).
% outputfile .. Base name for output files (string).
% 
% Usage:
% unwrapandnormalize(magfile,rawmagfile,phasefile,outputfile)

% Options
showplots=false;
subtractB0=false;
subtractmeanphase=false;

% IO
M=load_untouch_nii(magfile);
m = squeeze(M.img);

P=load_untouch_nii(phasefile); 
p=squeeze(P.img)/1000;

Mraw=load_untouch_nii(rawmagfile);

p=reshape(p,size(p,1),size(p,2),size(p,3),size(p,4)/4,4); % Unroll 5th dimension
p3d=reshape(p,size(p,1),size(p,2),[]); % Flatten into 3D

Miji; % Open ImageJ-Matlab interface
MIJ.createImage('phase3D',p3d,true)
MIJ.createImage('magnitude',m,true)

fprintf('\nRunning Phase Tools...')
fprintf('\nSet volume depth to %g',size(p,3))
fprintf('\nSet time steps to %g',size(p,4))
fprintf('\nSet images in parallel to %g',size(p,5))
fprintf(['\nWhen you are satisfied with the unwrap,' ...
         'return to Matlab and press any key to capture the image'])
fprintf('\nWaiting 5 seconds...\n')
pause(5);
MIJ.run('Phase Tools');
% TODO catch error here

pause(); % wait
unwrapped=MIJ.getCurrentImage(); % capture the active ImageJ image window
MIJ.exit; % close ImageJ

unwrapped=reshape(unwrapped,size(p)); % reshape into 5D
px=unwrapped(:,:,:,:,1); % x-sensitization
py=unwrapped(:,:,:,:,2); % y-sensitization
pz=unwrapped(:,:,:,:,3); % z-sensitization
pb=unwrapped(:,:,:,:,4); % B0-sensitization (magnetic field)

% normalize by subtracting mean phase across each slice
rx=px;
ry=py;
rz=pz;
rb=pb;
if subtractmeanphase
	for z=1:size(px,3),
		for t=1:size(px,4)
			% subtract background field from each slice
			nx(:,:,z,t)=px(:,:,z,t)-pb(:,:,z,t)*subtractB0;
			ny(:,:,z,t)=py(:,:,z,t)-pb(:,:,z,t)*subtractB0;
			nz(:,:,z,t)=pz(:,:,z,t)-pb(:,:,z,t)*subtractB0;
			nb(:,:,z,t)=pb(:,:,z,t);
			% find means of each field-corrected slice (ignoring background 0's)
			mx(z,t)=mean(nonzeros(nx(:,:,z,t)));
			my(z,t)=mean(nonzeros(ny(:,:,z,t)));
			mz(z,t)=mean(nonzeros(nz(:,:,z,t)));
			mb(z,t)=mean(nonzeros(nb(:,:,z,t)));
			% normalize by subtracting mean from non-background elements,
			% so that the new mean for every slice is 0
			rx(:,:,z,t)=nx(:,:,z,t)-mx(z,t)*(nx(:,:,z,t)~=0);
			ry(:,:,z,t)=ny(:,:,z,t)-my(z,t)*(ny(:,:,z,t)~=0);
			rz(:,:,z,t)=nz(:,:,z,t)-mz(z,t)*(nz(:,:,z,t)~=0);
			rb(:,:,z,t)=nb(:,:,z,t)-mb(z,t)*(nb(:,:,z,t)~=0);
			% find means of normalized slices (should all be 0)
			mnx(z,t)=mean(nonzeros(rx(:,:,z,t)));
			mny(z,t)=mean(nonzeros(ry(:,:,z,t)));
			mnz(z,t)=mean(nonzeros(rz(:,:,z,t)));
			mnb(z,t)=mean(nonzeros(rb(:,:,z,t)));
		end
	end
end

if showplots
	f1=figure; surfl(mx); title('X');
	f2=figure; surfl(my); title('Y');
	f3=figure; surfl(mz); title('Z');
	f4=figure; surfl(mb); title('B0');

	f5=figure; surfl(mnx); title('X-B0-mean');
	f6=figure; surfl(mny); title('Y-B0-mean');
	f7=figure; surfl(mnz); title('Z-B0-mean');
	f8=figure; surfl(mnb); title('B0-mean');
end

% Save NIFTI and Analyze (used in FSL, etc.)
% For this step, reassemble a combined magnitude/phase volume sequence
% from the raw magnitude data (4D) and unwrapped phase data (4D)
unwrappedFlt = reshape(...
					single(unwrapped),...
					size(P.img)...
					);
unwrappedInt=int32(unwrappedFlt*1000); % back to millirads in int32
nVolumes=size(unwrappedFlt,4);
combinedImg(:,:,:,1:nVolumes)=int32(Mraw.img);
combinedImg(:,:,:,(nVolumes+1):nVolumes*2)=unwrappedInt;

anaOut = make_ana(combinedImg,P.hdr.dime.pixdim(2:4));
niiOut = make_nii(combinedImg,P.hdr.dime.pixdim(2:4));
save_untouch_nii(anaOut,[outputfile '_unwrap.img']); % ANALYZE .img/.hdr
        save_nii(niiOut,[outputfile '_unwrap.nii']); % NIFTI .nii

% Save Philips PAR-REC (for Sinkus MRE software)
Analyze2PARREC([outputfile '_unwrap.img'])

% Save .MAT file (used in Matlab, Octave, MRE/Wave)
save([outputfile '_unwrap.mat'],'px','py','pz','rx','ry','rz','m')
save([outputfile '_workspace.mat'],'M','Mraw','P','unwrapped','niiOut')
