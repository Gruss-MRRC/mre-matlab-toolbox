function getUmap(im_mag,im_phase,sliceNum)

for j = 1:size(im_mag,1),
    for k = 1:size(im_mag,2),
        if (im_mag(j,k,4,1)>500),
            f1 = fftshift(fft(squeeze(im_phase(j,k,sliceNum,:,1))));
            f2 = fftshift(fft(squeeze(im_phase(j,k,sliceNum,:,2))));
            f3 = fftshift(fft(squeeze(im_phase(j,k,sliceNum,:,3))));
            f1_img(j,k) = abs(f1(6));
            f2_img(j,k) = abs(f2(6));
            f3_img(j,k) = abs(f3(6));
            ang1 = angle(f1);
            ang2 = angle(f2);
            ang3 = angle(f3);
            ang1_img(j,k) = ang1(6);
            ang2_img(j,k) = ang2(6);
            ang3_img(j,k) = ang3(6);
        else
            f1_img(j,k) = 0;
            f2_img(j,k) = 0;
            f3_img(j,k) = 0;
            ang1_img(j,k) = 0;
            ang2_img(j,k) = 0;
            ang3_img(j,k) = 0;
        end
        if (f2_img(j,k) > 3),
            %f2_img(j,k) = 0;
        end
        if (f3_img(j,k) > 3),
            %f3_img(j,k) = 0;
        end
    end
end

se = strel('disk',2);
BW1 = bwlabeln(imerode(im_mag(:,:,sliceNum,1)>50,se));
reg1 = regionprops(BW1,'Area');
[maxArea maxAreaIndex] = max([reg1.Area]);
BW2 = imdilate(BW1 == maxAreaIndex,se);

figure,subplot(1,3,1),imagesc(BW2.*f1_img)
subplot(1,3,2),imagesc(BW2.*f2_img)
subplot(1,3,3),imagesc(BW2.*f3_img)

figure,subplot(1,3,1),imagesc(BW2.*ang1_img)
subplot(1,3,2),imagesc(BW2.*ang2_img)
subplot(1,3,3),imagesc(BW2.*ang3_img)

f = bandpassfilter([size(f1_img,1),size(f1_img,2)], 0.03, 0.2, 5); %f(49,49)=1;
f1_img_filt = abs(ifft2(fftshift(fft2(f1_img)).*f));
f2_img_filt = abs(ifft2(fftshift(fft2(f2_img)).*f));
f3_img_filt = abs(ifft2(fftshift(fft2(f3_img)).*f));

curl1(1:size(im_mag,1),1:size(im_mag,2))=0;
curl2 = curl1;
curl3 = curl1;
for j = 2:size(im_mag,1)-1,
    for k = 2:size(im_mag,2)-1,
        curl1(j,k) = f1_img_filt(j-1,k) + f1_img_filt(j+1,k) + f1_img_filt(j,k-1) + f1_img_filt(j,k+1) - 4 * f1_img_filt(j,k);
        curl2(j,k) = f2_img_filt(j-1,k) + f2_img_filt(j+1,k) + f2_img_filt(j,k-1) + f2_img_filt(j,k+1) - 4 * f2_img_filt(j,k);
        curl3(j,k) = f3_img_filt(j-1,k) + f3_img_filt(j+1,k) + f3_img_filt(j,k-1) + f3_img_filt(j,k+1) - 4 * f3_img_filt(j,k);
    end
end

figure,subplot(1,3,1),imagesc(BW2.*curl1,[-.2 .2])
subplot(1,3,2),imagesc(BW2.*curl2,[-.2 .2])
subplot(1,3,3),imagesc(BW2.*curl3,[-.2 .2])
