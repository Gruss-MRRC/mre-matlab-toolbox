function frames = visualizeMotion(P,M,sliceN,exagFactor)
% Visualize 3D motion of a matrix slice
%
% Usage:
% frames = visualizeMotion(P,M,sliceN,exagFactor)
% 
% P ........... phase matrix; 5D image (xi,xj,xk,t,dxi)
% M ........... magnitude matrix, 3D image (xi,xj,xk). Ignore if [].
% sliceN....... slice number
% exagFactor .. scale motion by this amt
% frames ...... captured frames (play back with `movie(frames)`)

scale=exagFactor; %m/rad*exaggeration;
Sps=squeeze(P(:,:,sliceN,:,:))*scale;% scale
Sps((Sps==0))=NaN; 
width=size(Sps,1);
height=size(P,3);
N=size(Sps,3); % number of samples (time)
x=1:width;y=x;
[X,Y]=meshgrid(x,y);
Z=ones(width,width)*sliceN;

newsize = [floor(width/1.25) floor(width/1.25)];
X = imresize(X,newsize);
Y = imresize(Y,newsize);
Z = imresize(Z,newsize);

h = figure;
hold on;
colormap bone

if ~isempty(M)
	projectMagnitudeOntoWalls(h,M,false,sliceN)
end

axis([1 width 1 width -height*.25 height*1.5])
axis off
view([0.25 0.25 0.50])

for p=1:N;
    U = Sps(:,:,p,1); 
    V = Sps(:,:,p,2); 
    W = Sps(:,:,p,3); 
    U = imresize(U,newsize);
    V = imresize(V,newsize);
    W = imresize(W,newsize);
    %s = mesh(X+U,Y+V,Z+W);
   	s = surf(X+U,Y+V,Z+W);
    %set(s,'CData',sqrt(U.^2+V.^2+W.^2).^5);
    %set(s,'CData',W/2/pi);
    %set(s,'CDataMapping','scaled')
    %set(s,'AlphaData',.9*double(mag>100));
    set(s,'FaceColor',[.3 .3 .3]);
    set(s,'FaceAlpha',0.7);
    %set(s,'FaceAlpha','flat');
    %set(s,'EdgeColor','none');
    set(s,'EdgeColor',[0.6 0.8 1]);
  	F(p) = getframe(h);
  	delete(s)
end
frames = F;




