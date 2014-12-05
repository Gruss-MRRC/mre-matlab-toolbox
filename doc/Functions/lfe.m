function elastogram = lfe(img,mf,fov,scals,ws,norient,hm)
% Create elastogram from MRE phase data using local-frequency 
% estimate (LFE) algorithm
%
% Usage: 
% elastogram = lfe(img,mf,fov,scals,ws,norient,hm)
%
% img ......... 2 or 3D matrix, where the 3rd dim is time [radians]
% mf .......... mechanical frequency [Hz]
% fov ......... field-of-view (vector of length 2), [m]
% scals ....... filter spectrum max, log 2, e.g., spec = 2^[1 ... scals]
% ws .......... filter radius for smooting [px]
% norient ..... number of orientations
% hm .......... harmonic number
% 
% elastogram .. 2D matrix
%
% Authors:
% Ported from Java by Alex Krolick
% Original "MREJ" ImageJ plugin: 
% http://rsb.info.nih.gov/ij/plugins/mrej/index.html
% doi:10.1016/j.compbiomed.2013.04.005
%
% NOTE: Experimental! Output requires verification. As a direct port of the Java source code
% into Matlab, very little analysis has been done on the validity of this program.


  [rn, cn, sn] = size(img);
  %img=(img-min(img(:)))/(max(img(:))-min(img(:))); % normalize on range 0..1
  img = nthHarmonic(img,hm);
  imgfft=fftshift(fft2(img));
  switch scals<8
  	case true
  	 tb=2.^(1:scals);
 	   bw=2^(3/2);
 	   co=2^(1/2);
 	  case false
 	   tb=2.^((1:scals)/2);
     bw=2;
     co=2^(1/4);
  end
  thetasd = norient;
  rh=fov(1)/size(img,1);
  cl=fov(2)/size(img,2);
  mk=max(fov);
 	rho0=1/mk;
 	rho=rho0*tb;
 	slen=length(tb);
 	denom=zeros(size(img));
 	numer=denom;
 	for m=1:length(tb),
 	  qi{m}=zeros(size(imgfft));
 	  for n=1:norient,
 	    imgft=imgfft.*...
 	      lgnFilt(bw,[rn,cn],rho(m),pi/norient*n,norient,fov);
 	    imgft=fftshift(imgft);
 	    gt=ifft2(imgft);
 	    qi{m}=qi{m}+abs(gt);
 	  end
 	  denom=denom+qi{m};
 	  numer=numer+co*tb(m)*qi{m};
 	end
 	
 	lfeOut=rho0*numer./denom;
 	lfeOut=medfilt2(lfeOut,'indexed', [ws ws]);
 	lfeOut(lfeOut==0)=Inf;
 	elastogram=(mf./lfeOut).^2;
  
function lgn = lgnFilt(bw,filtsize,crho,theta,thetasigma,fov)
% Lognormal quadrature filter
% 
  cb = 4/log(2)/bw^2;
  rn = filtsize(1); cn = filtsize(2);
  xrange=linspace(-0.5*cn/fov(1),0.5*cn/fov(1),cn);
  yrange=linspace(-0.5*rn/fov(2),0.5*rn/fov(2),rn);
  [x y]=meshgrid(xrange,yrange);
  if thetasigma == 2
    thetasigma = pi / thetasigma / 2;
  elseif thetasigma == 4
    thetasigma = pi / thetasigma / 1.75;
  elseif thetasigma == 6
     thetasigma = pi / thetasigma / 1.5;
  end
  rd=zeros(rn,cn);
  lgn=ones(rn,cn);
  for ii=1:rn, 
    for jj=1:cn,
      if ii<length(x) && jj<length(y)
        rd(ii,jj)=sqrt(x(ii,jj)^2+y(ii,jj)^2);
      elseif ii>round(rn/2)-1 && ii<round(rn/2)+1 && jj>round(cn/2)+1 && jj<round(cn/2)+2,
        rd(ii,jj)=1;
        lgn(ii,jj)=0;
      end
    end
  end
  rd=log(rd/crho).^2;
  lgn=exp(-cb*rd).*lgn;
  T=atan2(-y,x); %rdtheta
  ds=sin(T)*cos(theta)-cos(T)*sin(theta);
  dc=cos(T)*cos(theta)+sin(T)*sin(theta);
  dtheta=abs(atan2(ds,dc));
  angularlgn=exp((-dtheta.^2)/(2*thetasigma^2));
  lgn=lgn.*angularlgn;
  
   
  
  
  
  
  
  
    
