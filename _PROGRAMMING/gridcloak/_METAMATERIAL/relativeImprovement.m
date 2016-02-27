function relimp = relativeImprovement(dphix,dphiy,alpha,beta,varargin)
%% plot optimized surface and raytracing results

%% input parsing
p = inputParser;

p.addParameter('n1',1,@isnumeric);              % starting ref. index
p.addParameter('n2',1.5,@isnumeric);            % end ref. index
p.addParameter('xlim',500,@isnumeric);          % x-size of unit cell
p.addParameter('ylim',50,@isnumeric);           % y-size of unit cell
p.addParameter('sx',0.9,@isnumeric);            % x-scaling factor
p.addParameter('sy',0.9,@isnumeric);            % y-scaling factor
p.addParameter('safetyFactor',0.9,@isnumeric);  % additional safety factor
p.addParameter('h0',100,@isnumeric);  
p.addParameter('xmin',-500,@isnumeric);  
p.addParameter('ymin',-50,@isnumeric);  
p.addParameter('numx',50,@isnumeric);  
p.addParameter('numy',5,@isnumeric);  

p.parse(varargin{:});

n1 = p.Results.n1;
n2 = p.Results.n2;
xlim = p.Results.xlim;
ylim = p.Results.ylim;
sx = p.Results.sx;
sy = p.Results.sy;
h0=p.Results.h0;
xmin=p.Results.xmin;
ymin=p.Results.ymin;
numx=p.Results.numx;
numy=p.Results.numy;

%% directions
v=[sind(beta) sind(alpha)*cosd(beta) cosd(alpha)*cosd(beta)];
n=[0 0 1];
%% orig
[x_orig,y_orig]=meshgrid(linspace(xmin+(xlim-xmin)/(2*numx),xlim-(xlim-xmin)/(2*numx),numx),linspace(ymin+(ylim-ymin)/(2*numy),ylim-(ylim-ymin)/(2*numy),numy));
%% tansmission
[v_trans]=Snell(dphix,dphiy,v,n,alpha,beta,'n1',n1,'n2',n2);
%% calculation of all the points
%points on Surface
H0=repmat(-h0,numy*numx,1);
SP=[x_orig(:),y_orig(:),H0(:)];
t=-H0./squeeze(v_trans(:,3));
T=repmat(t(:),1,3);
Schnittpunktebene=SP+v_trans.*T;

numOnContactGrid = length(find((1-sx)*xlim>mod((xlim+Schnittpunktebene(:,1)),(2*xlim)) |sx*xlim+xlim<mod((xlim+Schnittpunktebene(:,1)),(2*xlim))|(1-sy)*ylim>mod((ylim+Schnittpunktebene(:,2)),(2*ylim)) |sy*ylim+ylim<mod((ylim+Schnittpunktebene(:,2)),(2*ylim))));
np=size(Schnittpunktebene,1);
%% calculation of the relative Improvement
f=1-sx*sy;
N0=np;
N=np-numOnContactGrid;
relimp=N/(N0*(1-f))-1;
end