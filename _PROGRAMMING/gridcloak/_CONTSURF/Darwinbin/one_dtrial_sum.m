
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
xlim=500;
ylim=50;
nx=50;
ny=5;
h0=50;
%Parameter zum optimieren der Oberfläche
y0=40;
rx1=50;
rx2=500;
ry1=5;
ry2=50;
n=1.5;
xmax=500;
ymax=50;
c1=1.5;
c2=0.25;

[X,Y] = meshgrid(linspace(0,xlim,nx),linspace(0,ylim,ny));
Z0=-(c1*sqrt(y0*y0+rx1/(1-1/(n))*(2*(xmax-X)-((xmax-X).^(2))/(rx2)))+c2*sqrt(y0*y0+(ry1)/(1-1/(n))*(2*(ymax-Y)-((ymax-Y).^(2))/(ry2))))+( c1+c2)*y0-h0;
plotSurface(X,Y,Z0);