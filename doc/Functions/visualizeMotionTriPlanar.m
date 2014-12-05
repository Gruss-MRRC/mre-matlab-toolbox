function frames = visualizeMotionTriPlanar(P,sliceN,exagFactor)
% Visualize motion of a matrix slice in 3D
%
% Usage:
% frames = visualizeMotionTriPlanar(P,sliceN,exagFactor)
% 
% P ........... phase matrix; 5D image (xi,xj,xk,t,dxi)
% sliceN....... slice number
% exagFactor .. scale motion by this amt
% frames ...... captured frames (play back with `movie(frames)`)

scale=exagFactor; %m/rad*exaggeration;
Sps=squeeze(P(:,:,sliceN,:,:))*scale;
%Sps(isnan(Sps))=0;
Sps((Sps==0))=NaN; % force transparent background pixels
width=size(Sps,1);
height=size(P,3);
x=1:width;y=x;
[X,Y]=meshgrid(x,y);
Z=zeros(width,width);
N=size(Sps,3); % number of samples (time)

h = figure;
hold on;
colormap bone

% subtract DC components
subtractDC=false;
if subtractDC
for kk=1:3, dc{kk} = nthHarmonic(squeeze(Sps(:,:,:,kk)),1); end
for jj=1:3, for ii=1:N, 
	Sps(:,:,ii,jj) = Sps(:,:,ii,jj)-dc{jj};
end, end
end

ntnans=Sps(~isnan(Sps(:)));
cAx=[min(ntnans(:)) max(ntnans(:))];
%cAx=2*[-pi pi];

U = Sps(:,:,1,1); 
V = Sps(:,:,1,2); 
W = Sps(:,:,1,3); 

sp1 = subaxis(1,3,1, 'Spacing', 0.00, 'Padding', 0, 'Margin', 0);
s1 = surf(X,Y,Z);
set(s1,'CData',U);
set(s1,'CDataSource','U');
set(s1,'XDataSource','X');
set(s1,'FaceAlpha',1);
set(s1,'EdgeColor','none');
set(s1,'FaceColor','interp')
axis([1/7*width 6/7*width 1/7*width 6/7*width -1 1])
axis off
view(sp1,[0.05 0.10 -0.04])
caxis(cAx)
title('X Sensitization')

sp2 = subaxis(1,3,2, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
s2 = surf(X,Y,Z);
set(s2,'CData',V);
set(s2,'CDataSource','V');
set(s2,'YDataSource','Y');
set(s2,'FaceAlpha',1);
set(s2,'EdgeColor','none');
set(s2,'FaceColor','interp')
axis([1/7*width 6/7*width 1/7*width 6/7*width -1 1])
axis off
view(sp2,[0.05 0.10 -0.04])
caxis(cAx)
title('Y Sensitization')

sp3 = subaxis(1,3,3, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
s3 = surf(X,Y,Z);
set(s3,'CData',W);
set(s3,'CDataSource','W');
set(s3,'ZDataSource','Z');
set(s3,'FaceAlpha',1);
set(s3,'EdgeColor','none');
set(s3,'FaceColor','interp')
axis([1/7*width 6/7*width 1/7*width 6/7*width -1 1])
axis off
view(sp3,[0.05 0.10 -0.04])
caxis(cAx)
title('Z Sensitization')

for p=1:N;
    U = Sps(:,:,p,1); 
    V = Sps(:,:,p,2); 
    W = Sps(:,:,p,3); 
	refreshdata(h,'caller')
  	F(p) = getframe(h);
end
frames = F;






