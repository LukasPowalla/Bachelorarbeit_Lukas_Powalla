
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%% Parameter
h0=350;
numx=50;
numy=5;
numrayx=200;
numrayy=20;
xlim=500;
ylim=50;
xmin=-xlim;
ymin=-ylim;
n1=1;
n2=1.5;
safetyFactor=0.9;
sx=0.9;
sy=0.9;
orig_z=600;
%% Constructing of Surface
tr=constructSurface('h0',h0,'numx',numx,'numy',numy,'xlim',xlim,'ylim',ylim,'sx',sx,'sy',sy,'n1',n1,'n2',n2,'safetyFactor',safetyFactor);
%% annual intensity
nalpha=5;
nbeta=5;
[alphas,betas] = meshgrid(linspace(0,80,nalpha),linspace(0,80,nbeta));
P_values=zeros(nalpha,nbeta);
for i=1:nalpha
    for k=1:nalpha
    P_values(i,k)=relativeImprovement(tr,alphas(1,i),betas(k,1),'n1',n1,'n2',n2,'h0',h0,'numrayx',numrayx,'numrayy',numrayy,'xlim',xlim,'ylim',ylim,'sx',sx,'sy',sy,'orig_z',orig_z,'xmin',xmin,'ymin',ymin,'safetyFactor',safetyFactor)+1
    end
end

Q = annualImprovement(alphas,betas,P_values)

%% plotting
%[A,B]=[alphas,betas];
% surf(A,B,improvement);
% title(sprintf('max. improvement = %2.3f',max(improvement(:))));
% xlabel('alpha')
% ylabel('beta')
% zlabel('improvement');
% shading interp;
% colorbar;
