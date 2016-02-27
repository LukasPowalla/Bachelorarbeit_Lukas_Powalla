function [ fitness ] = getFitness( a,b,y0 )
%calulate fitness of sqrt(a*f(x)^2+b*g(y)^2) with uses ANALYTICAL APPROXIMATION 
%   Detailed explanation goes here


%% Parameter zum optimieren der Oberfläche
xlim=500;
ylim=50;
nx=50;
ny=5;
h0=60;
rx1=50;
rx2=500;
ry1=5;
ry2=50;
n=1.5;
xmax=500;
ymax=50;
%% fitness calculation
[X,Y] = meshgrid(linspace(0,xlim,nx),linspace(0,ylim,ny));
Z0=-sqrt(a*(y0*y0+rx1*(2*(xmax-X)-(xmax-X).^2/rx2)/(1-1/n))+b*(y0*y0+ry1*(2*(ymax-Y)-(ymax-Y).^2/ry2)/(1-1/n)))+sqrt(2*y0*y0)-h0;
fitness = plotSurface(X,Y,Z0);
end

