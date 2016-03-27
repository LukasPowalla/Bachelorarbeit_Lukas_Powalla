function [ Q ] = getQofFreeformwithParameter( varargin )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%% input parsing
p = inputParser;

p.addParameter('xlim',500,@isnumeric);              % starting ref. index
p.addParameter('ylim',50,@isnumeric);            % end ref. index
p.addParameter('nx',50,@isnumeric);          % x-size of unit cell
p.addParameter('ny',5,@isnumeric);           % y-size of unit cell
p.addParameter('rx1',50,@isnumeric);            % x-scaling factor
p.addParameter('rx2',500,@isnumeric);            % y-scaling factor
p.addParameter('ry1',5,@isnumeric);  % additional safety factor
p.addParameter('ry2',50,@isnumeric);
p.addParameter('a',2.07,@isnumeric);
p.addParameter('b',2.07,@isnumeric);
p.addParameter('y0x',50,@isnumeric);
p.addParameter('y0y',5,@isnumeric);

p.parse(varargin{:});

xlim = p.Results.xlim;
ylim = p.Results.ylim;
nx = p.Results.nx;
ny = p.Results.ny;
rx1 = p.Results.rx1;
rx2 = p.Results.rx2;
ry1 = p.Results.ry1;
ry2=p.Results.ry2;
a=p.Results.a;
b=p.Results.b;
y0x=p.Results.y0x;
y0y=p.Results.y0y;
%% Constructing of Surface
[X0,Y0] = meshgrid(linspace(0,xlim,nx),linspace(0,ylim,ny));
Z00=-sqrt(a*freeform(rx1,rx2,y0x,X0).^2+b*freeform(ry1,ry2,y0y,Y0).^2);
[X,Y] = meshgrid(linspace(-xlim,xlim,2*nx-1),linspace(-ylim,ylim,2*ny-1));
Z0 = [[rot90(Z00(2:end,2:end),2); fliplr(Z00(:,2:end))] [flipud(Z00(2:end,:));Z00]];
%% annual intensity
nalpha=10;
nbeta=10;
[alphas,betas] = meshgrid(linspace(0,80,nalpha),linspace(0,80,nbeta));
P_values=zeros(nalpha,nbeta);
for i=1:nalpha
    for k=1:nalpha
    P_values(i,k)=relativeImprovement(X,Y,Z0,alphas(1,i),betas(k,1))+1;
    end
end

Q = annualImprovement(alphas,betas,P_values)
end

