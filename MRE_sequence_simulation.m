
AT = 5; % min, total acq time
w1 = 55;  % vibration frequency, in Hz
T_1 = 900;
T_2 = 65;
T_2s = 30;

rfTime = 20; % ms, for 90, 180 and delay after sequence
ninetyTime = 14;
exTime = 2;
SENSEfact = 2.0;
nPE = 100;
nSlices = 8;
nPhases = 8;
nDirs = 4;

% Spin echo calculations
for j = 1:10, % vibration cycles
    if mod(j,2) == 0,
        cycles1(j)=2*j;
        ex_cycles = 2; % extra for shifting time
    else
        cycles1(j)=2*j; % extra for shifting time
        ex_cycles = 4; % extra for shifting time
    end
    for k = 1:10, % EPI acceleration
        acc1(k)=1+2*k;
        TE_SE(j,k) = 500*(cycles1(j)+ex_cycles)/w1+ (0.77 * acc1(k))/2 + rfTime;
        sliceTR(j,k) = TE_SE(j,k) + (0.77 * acc1(k))/2 + ninetyTime/2 + exTime;
        TR_SE(j,k)=sliceTR(j,k)*nSlices;
        Ns(j,k) = AT*60*1000/(sliceTR(j,k)* nSlices * nPhases * nDirs * nPE / SENSEfact / acc1(k));
        FA = acos(exp(-TR_SE(j,k)/T_1));
        SNR_SE(j,k) = sin(FA)*(1-exp(-TR_SE(j,k)/T_1))/(1-cos(FA)*exp(-TR_SE(j,k)/T_1))*exp(-TE_SE(j,k)/T_2)*(cycles1(j))*sqrt(Ns(j,k));
    end
end

BW = 145; %Hz
rfTime = 2;
acqTime = 7; %ms, for 145 Hz BW  (3 pixel WFS)
exTime = 1; % ms, for slice rephase

% Gradient echo calculations
for j = 1:5, % vibration cycles
    cycles1(j) = j;
    for k = 1:5, % vibration freq
        w2(k)=w1+30*(k-1);
        phase1 = pi * (1 - w1/w2(k));
        data1 = sin(2*pi*(1:1000)/1000).*sin(2*pi*(1:1000)/1000*w1/w2(k) + phase1);
        effic(k) = sum(data1)*w1/w2(k); % integrate over one period of the MEG
        TE(j,k) = rfTime/2 + 1000*cycles1(j)/w2(k) + acqTime/2 + exTime;
        TR(j,k) = TE(j,k) + rfTime/2 + acqTime/2 + exTime;
        Ernst(j) = acos(exp(-TR(j,k)/T_1));
        T1w(j) = (1-exp(-(TR(j,k))/T_1))*sin(Ernst(j))/(1-cos(Ernst(j))*exp(-(TR(j,k))/T_1));
        Ns(j,k) = AT*60*1000/(TR(j,k)* 1.3*nSlices * nPhases * nDirs * nPE / SENSEfact);
        SNR(j,k) = exp(-TE(j,k)/T_2s) * cycles1(j) * effic(k)/500 * sqrt(2200/BW)*sqrt(Ns(j,k)*nSlices)*T1w(j);
    end
end