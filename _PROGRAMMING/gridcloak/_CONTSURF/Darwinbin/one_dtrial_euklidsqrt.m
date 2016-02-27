% plotting of the solution sqrt(a*f(x)^2+b*g(y)^2) which uses ANALYTICAL APPROXIMATION 
%% allgemeine Parameter
xlim=500;
ylim=50;
nx=50;
ny=5;
%% Parameter zum optimieren der Oberfläche
y0=50;
rx1=50;
rx2=500;
ry1=5;
ry2=50;
n=1.5;
a=4;
b=3;
%% plotting
[X,Y] = meshgrid(linspace(0,xlim,nx),linspace(0,ylim,ny));
Z0=-sqrt(a*(y0*y0+rx1*(2*(xlim-X)-(xlim-X).^2/rx2)/(1-1/n))+b*(y0*y0+ry1*(2*(ylim-Y)-(ylim-Y).^2/ry2)/(1-1/n)))+sqrt(2*y0*y0);
plotSurface(X,Y,Z0,0,0);