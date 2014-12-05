function probePhaseData(imageDataM,imageDataP,imageDataS)

%b0valsPos = [1 8 18 29 43 54 65 76 87 98];
b0valsPos = [1 2];
but=1;
figure(2),imagesc(imageDataM(:,:,1)),colormap gray
while (but == 1),
    figure(2)
    [y1 x1 but] = ginput(1);
    x1=round(x1);
    y1=round(y1);
    figure(2),imagesc(imageDataM(:,:,1)),colormap gray
    figure(3),hold off, plot(unwrap(squeeze(imageDataM(x1,y1,:))))
    d1 = unwrap(squeeze(imageDataM(x1,y1,:))) - mean(unwrap(squeeze(imageDataM(x1,y1,:))));
    f1 = fftshift(fft(d1));
    posMid = floor(size(squeeze(imageDataM),3)/2) + 2;
    figure(4),hold off,plot(abs(f1(posMid-1:end)))
    bkgd = mean(abs(f1(end-5:end)));
    ratio1 = (sum(abs(f1(posMid+1:end))-bkgd)/(abs(f1(posMid))-bkgd));
    fourierPeak = max([abs(f1(posMid))]);
    text(2,0.1*fourierPeak,['Ratio1 = ' num2str(ratio1)]);
    
    %if imageDataP ~= [],
        figure(3),hold on,plot(squeeze(imageDataP(x1,y1,:)),'r-')
        d2 = unwrap(squeeze(imageDataP(x1,y1,:))) - mean(unwrap(squeeze(imageDataP(x1,y1,:))));
        f2 = fftshift(fft(d2));
        figure(4),hold on,plot(abs(f2(posMid-1:end)),'r-')
        ratio2 = (sum(abs(f2(posMid+1:end)))/abs(f2(posMid)));
        text(2,0.15*fourierPeak,['Ratio2 = ' num2str(ratio2)]);
        fourierPeak = max([abs(f1(posMid)) abs(f2(posMid))]);
    %end
    %if imageDataS ~= [],
        figure(3),hold on,plot(squeeze(imageDataS(x1,y1,:)),'g-')
        d3 = unwrap(squeeze(imageDataS(x1,y1,:))) - mean(unwrap(squeeze(imageDataS(x1,y1,:))));
        f3 = fftshift(fft(d3));
        figure(4),hold on,plot(abs(f3(posMid-1:end)),'g-')
        ratio3 = (sum(abs(f3(posMid+1:end)))/abs(f3(posMid)));
        text(2,0.2*fourierPeak,['Ratio3 = ' num2str(ratio3)]);
        fourierPeak = max([abs(f1(posMid)) abs(f2(posMid)) abs(f3(posMid))]);
   % end
    
end