function im1_p = spiralUnwrap(im1_m,im1_p,alpha)
% Unwrap a phase image along a spiral path

%alpha = 1.5;
% im1_m = im1(:,:,1:8);
% im1_p = im1(:,:,9:16);
im1_p = (im1_p - 2048)*pi/2048;

phantom1 = (imfill((im1_m(:,:,1)>200),'holes'));
for k = 1:size(im1_p,3),
    im1_p(:,:,k) = im1_p(:,:,k) .* int16(phantom1);
end

figure,imagesc(im1_p(:,:,8))
pCenter = round(ginput(1));
%pCenter = [125 118];
phantom2 = uint8(phantom1);
for j = 1:size(im1_p,3),
    insidePhantom = 1;
    stepPos = pCenter;
    numSteps = 0;
    stepForward = 1;
    while insidePhantom,
        numSteps = numSteps + 1;
        for stepDir = 0:1,
            for step = 1:numSteps,
                if stepDir,
                    corrSteps = 0;
                    if abs(im1_p(stepPos(2),stepPos(1)+stepForward,j)+pi)<0.0001,
                        % bad data point (original data was zero)
                        im1_p(stepPos(2),stepPos(1)+stepForward,j) = im1_p(stepPos(2),stepPos(1),j);
                    end
                    phaseStep = im1_p(stepPos(2),stepPos(1),j) - im1_p(stepPos(2),stepPos(1)+stepForward,j);
                    while abs(phaseStep) > alpha*pi && corrSteps < 4,
                        corrSteps = corrSteps + 1;
                        phantom2(stepPos(2),stepPos(1)) = 0;
                        if (phantom1(stepPos(2),stepPos(1))>0 && phantom1(stepPos(2),stepPos(1)+stepForward)>0),
                            stepSign = sign(phaseStep);
                            im1_p(stepPos(2),stepPos(1)+stepForward,j) = im1_p(stepPos(2),stepPos(1)+stepForward,j) + stepSign * 2 * pi;
                            phantom2(stepPos(2),stepPos(1)) = 2;
                        else
                            phaseStep = 0;
                        end
                        phaseStep = im1_p(stepPos(2),stepPos(1),j) - im1_p(stepPos(2),stepPos(1)+stepForward,j);
                    end
                    %figure(5),imagesc(phantom2)
                    if (mod(stepPos(2),5) == 0 && mod(stepPos(1),5) == 0),
                        figure(6),imagesc(im1_p(:,:,j))
                    end
                    stepPos(1) = stepPos(1) + stepForward;
                else
                    corrSteps = 0;
                    if abs(im1_p(stepPos(2)+stepForward,stepPos(1),j)+pi)<0.0001,
                        % bad data point (original data was zero)
                        im1_p(stepPos(2)+stepForward,stepPos(1),j) = im1_p(stepPos(2),stepPos(1),j);
                    end
                    phaseStep = im1_p(stepPos(2),stepPos(1),j) - im1_p(stepPos(2)+stepForward,stepPos(1),j);
                    while abs(phaseStep) > alpha*pi && corrSteps < 4,
                        corrSteps = corrSteps + 1;
                        phantom2(stepPos(2),stepPos(1)) = 0;
                        if (phantom1(stepPos(2),stepPos(1))>0 && phantom1(stepPos(2)+stepForward,stepPos(1))>0),
                            stepSign = sign(phaseStep);
                            im1_p(stepPos(2)+stepForward,stepPos(1),j) = im1_p(stepPos(2)+stepForward,stepPos(1),j) + stepSign * 2* pi;
                            phantom2(stepPos(2),stepPos(1)) = 4;
                        else
                            phaseStep = 0;
                        end
                        phaseStep = im1_p(stepPos(2),stepPos(1),j) - im1_p(stepPos(2)+stepForward,stepPos(1),j);
                    end
                    %figure(5),imagesc(phantom2)
                    %figure(6),imagesc(im1_p(:,:,j))
                    stepPos(2) = stepPos(2) + stepForward;
                end % if stepDir
                if phantom1(stepPos(2),stepPos(1))==0 && phantom1(stepPos(2),stepPos(1)+stepForward)==0 && phantom1(stepPos(2),stepPos(1)+2*stepForward)==0 && phantom1(stepPos(2),stepPos(1)+3*stepForward)==0,
                    insidePhantom = 0;
                end
                if phantom1(stepPos(2),stepPos(1))==0 && phantom1(stepPos(2)+stepForward,stepPos(1))==0 && phantom1(stepPos(2)+2*stepForward,stepPos(1))==0 && phantom1(stepPos(2)+3*stepForward,stepPos(1))==0,
                    insidePhantom = 0;
                end
            end % numSteps
        end % step direction
        stepForward = stepForward * -1;
    end % while insidePhantom
end % j loop
    