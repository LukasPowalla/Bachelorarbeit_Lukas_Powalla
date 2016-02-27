function [ fitness ] = getFitness_freeform( a,b,y0 )
%calulate fitness of sqrt(a*f(x)^2+b*g(y)^2) with uses NUMERICAL APPROXIMATION 
%   Detailed explanation goes here


%% Parameter zum optimieren der Oberfläche
xlim=500;
ylim=50;
nx=50;
ny=5;
h0=110;
rx1=50;
rx2=500;
ry1=5;
ry2=50;
n=1.5;
xmax=500;
ymax=50;
%% plotten
[X,Y] = meshgrid(linspace(0,xlim,nx),linspace(0,ylim,ny));
Z0=-sqrt(a*freeform(rx1,rx2,y0,X).^2+b*freeform(ry1,ry2,y0,Y).^2)+sqrt(2*y0*y0)-h0;
fitness = plotSurface(X,Y,Z0);
end

