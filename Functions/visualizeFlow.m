function frames = visualizeFlow(M,V,CM,C)
% Visualize out-of-plane fluid flow as a moving 3D surface
% Usage:
% frames = visualizeFlow(M,V,CM,C)
% M .............. Magnitude, 3D matrix (3rd dimension is time)
% V .............. Vascular velocity, 3D matrix
% MC ............. CSF magnitude image, 3D matrix
% C .............. CSF velocity, 3D matrix
% frames ......... Matlab frame object (struct of plot captures)
%                  Play back with movie(frames)
% NOTE: The particular preprocessing steps may need to be 
%       changed for a different set of images

% VASCULAR IMG PREPROCESS
  M = double(squeeze(M)); 
  V = -(double(squeeze(V))-2047); % center about 0
  V = V.*double(M>300); % mask
  % The middle values are right for the aliased section
  % flip V inside out and unwrap
  U = V;
  tmp=U(:,:,1);
  U(:,:,1) = U(:,:,8);
  U(:,:,8)=tmp;
  U = unwrap(U/2048*pi,[],3)/pi*2048;
  tmp=U(:,:,1);
  U(:,:,1) = U(:,:,8);
  U(:,:,8)=tmp;
  U(U>2200)=2200; % cap peaks
  U(U<-2200)=-2200;

% CSF IMG PREPROCESS
  CM = double(squeeze(CM));
  C = -(double(squeeze(C-2048))); % center about 0
  C = C.*double(CM>100); % mask
  C = C/4; % scale for plotting

% COMPUTE NET FLOWS
  T = size(M,3);
  for t = 1:T;
	  u = U(:,:,t); netu(t) = sum(sum(u(~isnan(u(:)))));
	  c = C(:,:,t); netc(t) = sum(sum(c(~isnan(c(:)))));
  end

% MAKE PLOTS
  h = figure();
  set(h,'color','black');
  
  % Show brain anatomy as a planar countour map of grey/white/csf.
  sub1 = axes('Position',[0.10 0.10 0.70 0.80]);
    hold on
    view([180-37.5 45])
    zoom(2.5)
    scale = 2048;
    axis([0 size(M,1) 0 size(M,2) -scale scale])
    axis('off');
    caxis([-scale,scale]);
    colormap gray
    m = M(:,:,1);
    peakM = max(m(:));
    cc = contour((m./peakM));
  
  % Fluid flow visualization
  sub2 = axes('Position',[0.05 0.85 0.95 0.10]);
    hold on
    plot([1:3*T],[netu netu netu],'r:','LineWidth',2)
    plot([1:3*T],[netc netc netc],'y:','LineWidth',2)
    axis tight
    axis off
    legend('Blood','CSF','Location','EO')
    legend boxoff
    %whitebg([0 0 0])
  for t_=[0 T 2*T]
  for t = 1:T
    u = U(:,:,t);
    c = C(:,:,t);
    c = medfilt2(c,[2,2]);
    u = medfilt2(u,[2,2]);
    cdatu = double(-(u<0)+(u>0));
    cdatc = double(-(c<0)+(c>0));
    u(imdilate(abs(cdatu),strel('disk',1))==0)=NaN;
    c(imdilate(abs(cdatc),strel('disk',1))==0)=NaN;

    % Color upward flow differently than downward flows
    % Blood is red/blue
    % CSF is yellow/cyan
    axes(sub1)
      s1 = surf(u.*(cdatu== 1),...%'CData',cdatu*scale,...
		       'FaceColor','red',...
		       'EdgeColor','red'...
		       );
      s2 = surf(u.*(cdatu==-1),...%'CData',cdatu*scale,...
		       'FaceColor','blue',...
		       'EdgeColor','blue'...
		       );
      s3 = surf(c.*(cdatc== 1),...%'CData',cdatc*scale/2,...
		       'FaceColor','yellow',...
		       'EdgeColor','yellow'...
		       );
      s4 = surf(c.*(cdatc==-1),...%'CData',cdatc*scale/2,...
		       'FaceColor','cyan',...
		       'EdgeColor','cyan'...
		       );
		     
    % Plot flow data on a line chart in the background
    axes(sub2)
      pt1 = plot(t+t_,netu(t),'r.','MarkerSize',20);
      pt2 = plot(t+t_,netc(t),'y.','MarkerSize',20);
	
    fprintf('\r%0.0d',t)
    frames(t+t_) = getframe(h);
    %pause(0.5)
    delete(s1),delete(s2),delete(s3),delete(s4),delete(pt1),delete(pt2)
  end
  end
  fprintf('\n')
