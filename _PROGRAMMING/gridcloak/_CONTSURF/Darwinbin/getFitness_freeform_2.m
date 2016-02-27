function [ fitness ] = getFitness_freeform_2( a,h0,y0 )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
%Parameter zum optimieren der Oberfläche
xlim=500;
ylim=50;
nx=50;
ny=5;
rx1=50;
rx2=500;
ry1=5;
ry2=50;

[X,Y] = meshgrid(linspace(0,xlim,nx),linspace(0,ylim,ny));
Z0=-sqrt(a*freeform(rx1,rx2,y0,X).^2+a*freeform(ry1,ry2,y0,Y).^2)+sqrt(2*y0*y0)-h0;
fitness = plotSurface(X,Y,Z0);
end

