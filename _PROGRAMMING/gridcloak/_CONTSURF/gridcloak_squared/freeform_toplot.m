% This script constructs ans plots the freeform for a quarter-cell with a given Parameter set
%% Parameter-set
xlim=500;
ylim=50;
nx=50;
ny=5;
h0=0;
rx1=50;
rx2=500;
ry1=5;
ry2=50;
n=1.5;
xmax=500;
ymax=50;
a=1;
b=1;
y0x=75;
y0y=7.5;
%% constructing surface
[X,Y] = meshgrid(linspace(0,xlim,nx),linspace(0,ylim,ny));
Z0=-sqrt(a*freeform(rx1,rx2,y0x,X).^2+b*freeform(ry1,ry2,y0y,Y).^2)-h0;
%% plotting surface
plotSurface(X,Y,Z0);