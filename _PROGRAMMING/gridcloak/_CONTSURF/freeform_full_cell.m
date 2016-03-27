% This script constructs the full cell with given Parameter set 
%% Parameterset
xlim=500;
ylim=50;
nx=50;
ny=5;
rx1=50;
rx2=500;
ry1=5;
ry2=50;
n=1.5;
a=1.9;
b=2.0;
y0x=75;
y0y=70;
%% Constructing of Surface
[X0,Y0] = meshgrid(linspace(0,xlim,nx),linspace(0,ylim,ny));
Z00=-sqrt(a*freeform(rx1,rx2,y0x,X0).^2+b*freeform(ry1,ry2,y0y,Y0).^2);
[X,Y] = meshgrid(linspace(-xlim,xlim,2*nx-1),linspace(-ylim,ylim,2*ny-1));
Z0 = [[rot90(Z00(2:end,2:end),2); fliplr(Z00(:,2:end))] [flipud(Z00(2:end,:));Z00]];
%% Plotting
% alpha=0;
% beta=0;
% plotSurface(X,Y,Z0,alpha,beta,'xlim',xlim,'ylim',ylim);
%display(relativeImprovement(X,Y,Z0,alpha1,beta1,'xlim',xlim,'ylim',ylim));
%% annual intensity
nalpha=50;
nbeta=50;
[alphas,betas] = meshgrid(linspace(0,80,nalpha),linspace(0,80,nbeta));
P_values=zeros(nalpha,nbeta);
for i=1:nalpha
    for k=1:nalpha
    P_values(i,k)=relativeImprovement(X,Y,Z0,alphas(1,i),betas(k,1))+1;
    end
end

Q = annualImprovement(alphas,betas,P_values)