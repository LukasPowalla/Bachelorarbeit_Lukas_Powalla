function [ Q ] = getQwithParameter( varargin )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%% input parsing
p = inputParser;

p = inputParser;

p.addParameter('n1',1,@isnumeric);              % starting ref. index
p.addParameter('n2',1.5,@isnumeric);            % end ref. index
p.addParameter('xlim',500,@isnumeric);          % x-size of unit cell
p.addParameter('ylim',50,@isnumeric);           % y-size of unit cell
p.addParameter('xmin',-500,@isnumeric);  
p.addParameter('ymin',-50,@isnumeric);  
p.addParameter('sx',0.9,@isnumeric);            % x-scaling factor
p.addParameter('sy',0.9,@isnumeric);            % y-scaling factor
p.addParameter('safetyFactor',0.9,@isnumeric);  % additional safety factor
p.addParameter('h0',30,@isnumeric);   
p.addParameter('numx',50,@isnumeric);  
p.addParameter('numy',5,@isnumeric);  


p.parse(varargin{:});

n1 = p.Results.n1;
n2 = p.Results.n2;
xlim = p.Results.xlim;
ylim = p.Results.ylim;
sx = p.Results.sx;
sy = p.Results.sy;
safetyFactor=p.Results.safetyFactor;
h0=p.Results.h0;
numx=p.Results.numx;
numy=p.Results.numy;
xmin=p.Results.xmin;
ymin=p.Results.ymin;
%% surface construction
[dphix,dphiy]=constructSurface('n1',n1,'n2',n2,'h0',h0,'xlim',xlim,'ylim',ylim,'sx',sx,'sy',sy,'safetyFactor',safetyFactor,'numx',numx,'numy',numy);
%% annual intensity
nalpha=100;
nbeta=100;
[alphas,betas] = meshgrid(linspace(0,80,nalpha),linspace(0,80,nbeta));
P_values=zeros(nalpha,nbeta);
for i=1:nalpha
    for k=1:nalpha
    P_values(i,k)=relativeImprovement(dphix,dphiy,alphas(1,i),betas(k,1),'n1',n1,'n2',n2,'h0',h0,'xlim',xlim,'ylim',ylim,'sx',sx,'sy',sy,'xmin',xmin,'ymin',ymin,'safetyFactor',safetyFactor)+1;
    end
end

Q = annualImprovement(alphas,betas,P_values)
end

