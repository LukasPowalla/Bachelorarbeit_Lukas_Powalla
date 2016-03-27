% This script constructs the full cell with given Parameter set 
%% Parameter-set
xlim=50;
ylim=50;
nx=15;
ny=15;
rx1=5;
rx2=50;
ry1=5;
ry2=50;
n=1.5;
a=0.6;
b=1.7;
y0x=5.5;
y0y=6.0;
%% Constructing of Surface
[X0,Y0] = meshgrid(linspace(0,xlim,nx),linspace(0,ylim,ny));
Z00=-sqrt(a*freeform(rx1,rx2,y0x,X0).^2+b*freeform(ry1,ry2,y0y,Y0).^2);
[X,Y] = meshgrid(linspace(-xlim,xlim,2*nx-1),linspace(-ylim,ylim,2*ny-1));
Z0 = [[rot90(Z00(2:end,2:end),2); fliplr(Z00(:,2:end))] [flipud(Z00(2:end,:));Z00]];
%% Plotting
 alpha=0;
 beta=0;
 plotSurface(X,Y,Z0,alpha,beta,'xlim',xlim,'ylim',ylim);
% display(relativeImprovement(X,Y,Z0,alpha1,beta1,'xlim',xlim,'ylim',ylim));
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