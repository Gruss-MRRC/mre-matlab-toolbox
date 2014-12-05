function varargout = splitDimension(matrix,varargin)
% Split matrix along dimension
% If no dimension specified, assumes highest dimension
% Usage:
% [varargout] = splitDimension(matrix,<optional dimension number>)

% Example:
% >> matIn = ones(5,5,3,2)
% >> [mat1,mat2] = splitDimension(matIn);
% >> size(mat1), size(mat2)
%    5    5    3
%    5    5    3
% >> [mat1,mat2,mat3] = splitDimension(matIn,3);
% >> size(mat1), size(mat2), size(mat3)
%    5    5    2
%    5    5    2
%    5    5    2

warning('Warning: Unfinished function. SPLITDIMENSION only works on 5D matrices!') 
switch nargin
  case 1
    dim = length(size(matrix));
  case 2
    dim = varargin{1};
  otherwise
    errmsg = [nargin ' arguments given, expected 1 or 2'];
    error(errmsg)
end

for index = 1:size(matrix,dim)
  varargout{index} = matrix(:,:,:,:,index);
end

