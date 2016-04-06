%script to construct and raytrace a metamaterial
%% an. in deg. 
alpha1=0; %angles for optimisation
beta1=0;
alpha=0; %angles for ray tracing
beta=0;
%% parameter
n1 = 1;
n2 = 1.5;
xlim = 500;
ylim = 50;
sx = 0.9;
sy = 0.9;
safetyFactor = 0.9;
numx=50;
numy=10;
h0=100;
k0=1;
orig_z=200;
xmin=-xlim;
ymin=-ylim;
%% surface
[dphix,dphiy,x,y,z,dxz]=constructSurface('numx',numx,'numy',numy,'sx',sx,'sy',sy,'n1',n1,'n2',n2);
n=[0 0 1];
%% directions
v=[sind(beta) sind(alpha)*cosd(beta) cosd(alpha)*cosd(beta)];
dir=repmat(v,(2*numx-1)*(2*numy-1),1);
%% orig
x_orig=x;
y_orig=y;
z_orig=z;
orig=[x_orig(:) y_orig(:) z_orig(:)];
orig=orig-dir.*200;
%% tansmission
[v_trans]=Snell(dphix,dphiy,dxz,alpha,beta,'n1',n1,'n2',n2);
%% calculation of all the points
%points on Surface
SP=[x(:),y(:),z(:)];
t=-z_orig(:)./squeeze(v_trans(:,3));
T=repmat(t(:),1,3);
SP_ebene=SP+v_trans.*T;
%% PlotSurface

plotSurface(orig,SP,SP_ebene,x,y,z);