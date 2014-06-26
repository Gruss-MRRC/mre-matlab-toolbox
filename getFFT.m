function getFFT(theImg)

but=1;
while but == 1,
figure(2),imagesc(mean(theImg,3))
[y1 x1 but] = (ginput(1));y1=round(y1);x1=round(x1);
data1 = double(squeeze(mean(mean(theImg(x1-1:x1+1,y1-1:y1+1,:),1),2))-mean(squeeze(mean(mean(theImg(x1-1:x1+1,y1-1:y1+1,:),1),2))));
f1 = fftshift(fft(data1));
figure(3),plot(data1(1:30))
figure(4),plot(abs(f1(1:end)),'x-')
end