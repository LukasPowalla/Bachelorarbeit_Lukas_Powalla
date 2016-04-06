function [ dphix,dphiy,x,y,z,dxz ] = constructSurface( varargin )
%constructs Surface, especially the matrices dphix and dphiy
%% input parsing
p = inputParser;

p.addParameter('n1',1,@isnumeric);              % starting ref. index
p.addParameter('n2',1.5,@isnumeric);            % end ref. index
p.addParameter('xlim',500,@isnumeric);          % x-size of un1t cell
p.addParameter('ylim',50,@isnumeric);           % y-size of un1t cell
p.addParameter('numx',50,@isnumeric);
p.addParameter('numy',10,@isnumeric);
p.addParameter('sx',0.9,@isnumeric);            % x-scaling factor
p.addParameter('sy',0.9,@isnumeric);            % y-scaling factor
p.addParameter('safetyFactor',0.9,@isnumeric);  % additional safety factor
p.addParameter('k0',1,@isnumeric);


p.parse(varargin{:});


n1 = p.Results.n1;
n2 = p.Results.n2;
xlim = p.Results.xlim;
ylim = p.Results.ylim;
sx = p.Results.sx;
sy = p.Results.sy;
safetyFactor = p.Results.safetyFactor;
numx=p.Results.numx;
numy=p.Results.numy;
k0=p.Results.k0;
%% Parameterset
rx1=(1-sx*safetyFactor)*xlim;
rx2=xlim;
a=1;
y0x=rx1;
%% safetyfactor
sy=sy*safetyFactor;
%% Constructing of Surface
[X0,Y0] = meshgrid(linspace(0,xlim,numx),linspace(0,ylim,numy));
Z00=-a*freeform(rx1,rx2,y0x,X0);
[x,y] = meshgrid(linspace(-xlim,xlim,2*numx-1),linspace(-ylim,ylim,2*numy-1));
z = [[rot90(Z00(2:end,2:end),2); fliplr(Z00(:,2:end))] [flipud(Z00(2:end,:));Z00]];
%% calculation of the derivative dxz
dxz=gradient(z(1,:))./gradient(x(1,:));
dxz= repmat(dxz(1,:),[size(z(:,1),1) 1]);
dphiy=-k0.*n2.*sqrt(1-n1.^2./n2.^2.*sind(acosd(((dxz.^2+1).^(-1./2))).*sign(dxz)).^2).*sign(dxz).*(sy-1).*y.*((dxz.^2+1).^(-1./2)).*(abs(dxz).*dxz+sign(dxz)).*((1+sign(dxz).^2.*(sy-1).^2.*y.^2./(dxz.^2+1)./(abs(dxz).*x.*sx+sign(dxz).*z-abs(dxz).*x).^2.*(abs(dxz).*dxz+sign(dxz)).^2).^(-1./2))./(abs(dxz).*(sx-1).*x+sign(dxz).*z);
dphix=0;
end

