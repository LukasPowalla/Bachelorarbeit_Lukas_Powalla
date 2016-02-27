function [ output_args ] = plotSurface(tr1,orig,Schnittpunkte1,Schnittpunktebene,varargin)
%plots the surface
%   Detailed explanation goes here

%% input parsing
p = inputParser;

p.addOptional('alpha',0,@isnumeric);
p.addOptional('beta',0,@isnumeric);
p.addParameter('n1',1,@isnumeric);              % starting ref. index
p.addParameter('n2',1.5,@isnumeric);            % end ref. index
p.addParameter('xlim',500,@isnumeric);          % x-size of unit cell
p.addParameter('ylim',50,@isnumeric);           % y-size of unit cell
p.addParameter('sx',0.9,@isnumeric);            % x-scaling factor
p.addParameter('sy',0.9,@isnumeric);            % y-scaling factor
p.addParameter('safetyFactor',0.9,@isnumeric);  % additional safety factor
p.addParameter('h0',25,@isnumeric);
p.addParameter('numx',20,@isnumeric);
p.addParameter('numy',5,@isnumeric);

p.parse(varargin{:});

alpha = p.Results.alpha;
beta = p.Results.beta;
n1 = p.Results.n1;
n2 = p.Results.n2;
xlim = p.Results.xlim;
ylim = p.Results.ylim;
sx = p.Results.sx;
sy = p.Results.sy;
safetyFactor = p.Results.safetyFactor;
h0=p.Results.h0;
numx=p.Results.numx;
numy=p.Results.numy;

%% plot
np=size(orig,1);
figure(133)
trisurf(tr1,'FaceColor','interp');
hold on;
patch([xlim xlim -xlim -xlim sx*xlim sx*xlim xlim],[-ylim ylim ylim sy*ylim sy*ylim -ylim -ylim],zeros(1,7),'FaceColor','k','FaceAlpha',.5);
patch([-xlim -xlim xlim xlim -sx*xlim -sx*xlim -xlim],[ylim -ylim -ylim -sy*ylim -sy*ylim ylim ylim],zeros(1,7),'FaceColor','k','FaceAlpha',.5);
line([orig(:,1)';Schnittpunkte1(:,1)';Schnittpunktebene(:,1)'],[orig(:,2)';Schnittpunkte1(:,2)';Schnittpunktebene(:,2)'],[orig(:,3)';Schnittpunkte1(:,3)';Schnittpunktebene(:,3)'],'Color','k','LineWidth',0.3);
%line([orig(:,1)';orig(:,1)'],[orig(:,2)';orig(:,2)'],[orig(:,3)';zeros(size(orig(:,1)))'],'Color','k','LineWidth',0.3);
scatter3(Schnittpunktebene(:,1)',Schnittpunktebene(:,2)',zeros(1,np),'*r');
scatter3(sx*Schnittpunkte1(:,1)',sy*Schnittpunkte1(:,2)',zeros(1,np),'.g');
line([Schnittpunktebene(:,1)';sx*orig(:,1)'], [Schnittpunktebene(:,2)';sy*orig(:,2)'], zeros(2,np),'Color','r');
title('fresnel construction');
xlabel('x');
ylabel('y');
zlabel('z');
%axis equal;
hold off

end

