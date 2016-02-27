function [A,B] = cloneToEqualSize(A0,B0,varargin)
%% Clones matrices to matching size along dimension dim.
% By default, dim is first dimension with non-equal size.
% Output matrices have the same size along dimension dim. In case the
% sizes of A0 and B0 along dim are not divisble, the resulting cloned
% matrix is cropped to fit the size of the larger matrix of A0 and B0.
%
% Example:
%  A0 = [1 2; 3 4];  B0 = [1 2; 3 4; 5 6];
%  -> dim = 1 (because A0's and B0's row count differs)
%  -> A = [1 2; 3 4; 1 2];
%  -> B = [1 2; 3 4; 5 6];

p = inputParser;
p.addOptional('dim',0,@isnumeric);
p.parse(varargin{:});
dim = p.Results.dim;

if ~isscalar(dim) || ~isfinite(dim) || uint8(dim)~=dim || isnan(dim) || dim<0 || dim>3
    error('dim must be a positive integer scalar between 1 and 3');
end

sA = size(A0); 
sB = size(B0);

if length(sA)~=length(sB)
    error('A0 and B0 must have same dimensionality');
end

if dim==0
    dim = find(sA~=sB,1,'first');
    if isempty(dim)
        A = A0;
        B = B0;
        return;
    end
end

sA(dim) = 1;
sB(dim) = 1;
if any(sA~=sB)
    warning('sizes of A0 and B0 along remaining dimensions are not equal');
end

sA0 = size(A0,dim);
sB0 = size(B0,dim);

A = A0;
B = B0;
if sA0>sB0
    repvec = ones(size(sB));
    repvec(dim) = ceil(sA0/sB0);
    B = repmat(B0,repvec);
    C={};
    for m = 1:length(sB)
        s = sB(m);
        C{m} = 1:s;
    end
    C{dim} = 1:size(A0,dim);
    B = B(C{:});
elseif sB0>sA0
    repvec = ones(size(sA));
    repvec(dim) = ceil(sB0/sA0);
    A = repmat(A0,repvec);
    C={};
    for m = 1:length(sA)
        s = sA(m);
        C{m} = 1:s;
    end
    C{dim} = 1:size(B0,dim);
    A = A(C{:});
end