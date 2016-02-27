
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
xlim=500;
ylim=50;
nx=50;
ny=5;
h0=50;
%Parameter zum optimieren der Oberfläche
ampl = 200;
b=1;
c=0.5;

[X,Y] = meshgrid(linspace(0,xlim,nx),linspace(0,ylim,ny));
Z0 =(ampl*(c*X.^2.+b*Y.^2)/max(c*X(:).^2+b*Y(:).^2)-h0-ampl);

plotSurface(X,Y,Z0);

