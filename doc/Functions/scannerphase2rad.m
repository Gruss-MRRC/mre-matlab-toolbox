function radianphaseimage = scannerphase2rad(rawphaseimage)
% Convert units of phase image from MRI scanner (range 0...2^n) to radians
% (range -pi...pi)
% Takes: rawphaseimage (any dimension array, any numeric type)
% Returns: radianphaseimage (same dimension, type double float)
backgroundpixels=(rawphaseimage==0);
radianphaseimage = double(rawphaseimage)*pi/2048-pi;
radianphaseimage(backgroundpixels)=0;