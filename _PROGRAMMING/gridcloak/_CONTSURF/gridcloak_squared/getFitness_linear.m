function [ fitness ] = getFitness_linear( c1,c2,y0 )
%calulate fitness of a*f(x)+b*g(y) with uses ANALYTICAL APPROXIMATION 
%   Detailed explanation goes here


%% Parameter zum optimieren der Oberfläche
xlim=500;
ylim=50;
nx=50;
ny=5;
rx1=50;
rx2=500;
ry1=5;
ry2=50;
n=1.5;
xmax=500;
ymax=50;
%% plotting
[X,Y] = meshgrid(linspace(0,xlim,nx),linspace(0,ylim,ny));
Z0=-(c1*sqrt(y0*y0+rx1/(1-1/(n))*(2*(xmax-X)-((xmax-X).^(2))/(rx2)))+c2*sqrt(y0*y0+(ry1)/(1-1/(n))*(2*(ymax-Y)-((ymax-Y).^(2))/(ry2))));
fitness = plotSurface(X,Y,Z0);
end

