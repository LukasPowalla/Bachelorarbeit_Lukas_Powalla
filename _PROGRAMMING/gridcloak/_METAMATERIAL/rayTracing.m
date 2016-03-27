%script to construct and raytrace a metamaterial
%% an. in deg. 
alpha=20;
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
h0=30;
k0=1;
orig_z=200;
xmin=-xlim;
ymin=-ylim;
%% directions
v=[sind(beta) sind(alpha)*cosd(beta) cosd(alpha)*cosd(beta)];
dir=repmat(v,numx*numy,1);
%% orig
[x_orig,y_orig]=meshgrid(linspace(xmin+(xlim-xmin)/(2*numx),xlim-(xlim-xmin)/(2*numx),numx),linspace(ymin+(ylim-ymin)/(2*numy),ylim-(ylim-ymin)/(2*numy),numy));
z_orig=zeros((numy)*(numx),1);
orig=[x_orig(:) y_orig(:) z_orig];
orig=(-(orig_z-h0)/(cosd(alpha)*cosd(beta)))*dir+orig;
orig(:,3)=orig(:,3)-h0;

%% surface
[dphix,dphiy]=constructSurface('numx',numx,'numy',numy,'sx',sx,'sy',sy,'h0',h0,'n1',n1,'n2',n2);
n=[0 0 1];
%% tansmission
[v_trans]=Snell(dphix,dphiy,v,n,alpha,beta,'n1',n1,'n2',n2);
%% calculation of all the points
%points on Surface
H0=repmat(-h0,numy*numx,1);
SP=[x_orig(:),y_orig(:),H0(:)];
t=-H0./squeeze(v_trans(:,3));
T=repmat(t(:),1,3);
SP_ebene=SP+v_trans.*T;
%% PlotSurface

plotSurface(orig,SP,SP_ebene,'h0',h0);