
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
xlim=500;
ylim=50;
nx=50;
ny=5;
h0=50;
%Parameter zum optimieren der Oberfläche

b=500;
c=50;
d=1;
e=3;


[X,Y] = meshgrid(linspace(0,xlim,nx),linspace(0,ylim,ny));
Z0=-ampl*cos(X*pi/(2*d*b)).*cos(Y*pi/(2*c*e))-h0
plotSurface(X,Y,Z0);
