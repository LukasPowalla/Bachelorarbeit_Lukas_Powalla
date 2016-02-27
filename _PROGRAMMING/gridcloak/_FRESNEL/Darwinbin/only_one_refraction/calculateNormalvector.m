function [nx,ny,nz] = calculateNormalvector(x,y,varargin)
%calculates Normalvector for given meshgrid (x|y) with alpha=0,beta=0
%   for normal incidence
%% input parsing
p = inputParser;

p.addParameter('n1',1,@isnumeric);              % starting ref. index
p.addParameter('n2',1.5,@isnumeric);            % end ref. index
p.addParameter('h0',25,@isnumeric);             % minimum surface height
p.addParameter('sx',0.9,@isnumeric);            % x-scaling factor
p.addParameter('sy',0.9,@isnumeric);            % y-scaling factor

p.parse(varargin{:});


n1 = p.Results.n1;
n2 = p.Results.n2;
h0 = p.Results.h0;
sx = p.Results.sx;
sy = p.Results.sy;

%% calculaion of normalvector with meshgrid
h0=-h0;
nx=-(n2*(sx-1)*x)./sqrt(h0^2*(n1^2+n2^2)+2*h0*n1*n2*sqrt(h0^2+(sx-1)^2*x.^2+(sy-1)^2*y.^2)+(n1^2+n2^2)*((sx-1)^2*x.^2+(sy-1)^2*y.^2));
ny=-(n2*(sy-1)*y)./sqrt(h0^2*(n1^2+n2^2)+2*h0*n1*n2*sqrt(h0^2+(sx-1)^2*x.^2+(sy-1)^2*y.^2)+(n1^2+n2^2)*((sx-1)^2*x.^2+(sy-1)^2*y.^2));
nz=(n1*sqrt(h0^2+(sx-1)^2*x.^2+(sy-1)^2*y.^2)+h0*n2)./sqrt(h0^2*(n1^2+n2^2)+2*h0*n1*n2*sqrt(h0^2+(sx-1)^2*x.^2+(sy-1)^2*y.^2)+(n1^2+n2^2)*((sx-1)^2*x.^2+(sy-1)^2*y.^2));
nz=-nz;
end

