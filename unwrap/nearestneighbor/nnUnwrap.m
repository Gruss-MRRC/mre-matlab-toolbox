function [im1_p,centPos] = nnUnwrap(im1_m,im1_p,alpha,centPos,showFig,P_axes)
% Nearest-Neighbor Unwrap

%alpha = 1.5;
im1_p = double(im1_p - 2048)*pi/2048;

phantom1 = double(imfill((im1_m(:,:,1)>200),'holes'));
for k = 1:size(im1_p,3),
    im1_p(:,:,k) = im1_p(:,:,k) .* phantom1;
end

if isempty(centPos)
    imagesc(im1_p(:,:,1),'Parent',P_axes) % get center
    centPos = round(ginput(1));
end

%centPos = [125 118];
phantom2 = uint8(phantom1);
s1 = size(im1_p);

phaseCent = mean(mean(im1_p(centPos(2)-2:centPos(2)+2,centPos(1)-2:centPos(1)+2,1),1),2);
for j = 1:s1(3),
    stdPhase = std(nonzeros(im1_p(centPos(2)-1:centPos(2)+1,centPos(1)-1:centPos(1)+1,j)));
    if stdPhase > 1,
        % reacquire center, too close to phase alias edge
        imagesc(im1_p(:,:,j),'Parent',P_axes) % get center
        centPos = round(ginput(1));
    end
    insidePhantom = 0;
    stepPos = centPos;
    numSteps = 0;
    stepForward = 1;
    while insidePhantom < 500 && stepPos(2) < s1(1)-1 && stepPos(2) > 1 && stepPos(1) < s1(2)-1 && stepPos(1) > 1,
        numSteps = numSteps + 1;
        for stepDir = 0:1,
            for step = 1:numSteps,
                corrSteps = 0;
                if abs(im1_p(stepPos(2),stepPos(1)+stepForward,j)+pi)<0.0001,
                    % bad data point (original data was zero)
                    im1_p(stepPos(2),stepPos(1)+stepForward,j) = im1_p(stepPos(2),stepPos(1),j);
                end
                for jj = -1:2:1,
                    for kk = -1:2:1,
                        phaseStep = im1_p(stepPos(2),stepPos(1),j) - im1_p(stepPos(2)+jj,stepPos(1)+kk,j);
                        if (phantom1(stepPos(2),stepPos(1))>0 && phantom1(stepPos(2)+jj,stepPos(1)+kk)>0),
                            insidePhantom = 0;
                            while abs(phaseStep) > alpha*pi && corrSteps < 4,
                                corrSteps = corrSteps + 1;
                                phantom2(stepPos(2),stepPos(1)) = 0;
                                stepSign = sign(phaseStep);
                                im1_p(stepPos(2)+jj,stepPos(1)+kk,j) = im1_p(stepPos(2)+jj,stepPos(1)+kk,j) + stepSign * 2 * pi;
                                phantom2(stepPos(2),stepPos(1)) = 2;    
                                phaseStep = im1_p(stepPos(2),stepPos(1),j) - im1_p(stepPos(2)+jj,stepPos(1)+kk,j);
                            end % do phase correction
                        end % phantom1 > 0
                    end % kk
                end % jj
                %figure(5),imagesc(phantom2)
                if (showFig && mod(stepPos(2),5) == 0 && mod(stepPos(1),5) == 0),
                   figure(6),imagesc(im1_p(:,:,j))
                end
                
                if stepDir,
                   stepPos(1) = stepPos(1) + stepForward;
                else    
                    stepPos(2) = stepPos(2) + stepForward;
                end % if stepDir
                
                if phantom1(stepPos(2),stepPos(1))==0,
                    insidePhantom = insidePhantom + 1;
                end

            end % numSteps
        end % step direction
        stepForward = stepForward * -1;
    end % while insidePhantom
    
    phaseCentCurr = mean(mean(im1_p(centPos(2)-2:centPos(2)+2,centPos(1)-2:centPos(1)+2,j),1),2);
    
    if abs(phaseCentCurr - phaseCent) > 2.0,
        sgnPhase = sign(phaseCentCurr - phaseCent);
        im1_p(:,:,j) = im1_p(:,:,j) - (im1_p(:,:,j) ~= 0) * sgnPhase * 2 * pi;
    end
    phaseCent = mean(mean(im1_p(centPos(2)-2:centPos(2)+2,centPos(1)-2:centPos(1)+2,j),1),2);
end % j loop
    