% Creates bandpass Butterworth filter in two dimensions.
%
% BANDPASSFILTER - Constructs a band-pass butterworth filter
%
% usage: f = bandpassfilter(sze, cutin, cutoff, n)
% 
% where: sze    is a two element vector specifying the size of filter 
%               to construct.
%        cutin and cutoff are the frequencies defining the band pass 0 - 0.5
%        n      is the order of the filter, the higher n is the sharper
%               the transition is. (n must be an integer >= 1).
%
% The frequency origin of the returned filter is at the corners.
%
% See also: LOWPASSFILTER, HIGHPASSFILTER, HIGHBOOSTFILTER
%

% Peter Kovesi   pk@cs.uwa.edu.au
% Department of Computer Science & Software Engineering
% The University of Western Australia
%
% October 1999

function f = bandpassfilter(sze, cutin, cutoff, n)
    
    if cutin < 0 || cutin > 0.5 || cutoff < 0 || cutoff > 0.5
	error('frequencies must be between 0 and 0.5');
    end
    
    if rem(n,1) ~= 0 || n < 1
	error('n must be an integer >= 1');
    end
    
f = lowpassfilter(sze, cutoff, n) - lowpassfilter(sze, cutin, n);