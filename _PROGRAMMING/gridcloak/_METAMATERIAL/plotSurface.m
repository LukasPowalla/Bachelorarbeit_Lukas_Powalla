function [ ] = plotSurface(orig,Schnittpunkt,Schnittpunktebene,varargin)
%plots the surface
%  plots the surface :)

%% input parsing
p = inputParser;

p.addParameter('xlim',500,@isnumeric);          % x-size of unit cell
p.addParameter('ylim',50,@isnumeric);           % y-size of unit cell
p.addParameter('sx',0.9,@isnumeric);            % x-scaling factor
p.addParameter('sy',0.9,@isnumeric);            % y-scaling factor
p.addParameter('safetyFactor',0.9,@isnumeric);  % additional safety factor
p.addParameter('h0',25,@isnumeric);
p.addParameter('numx',20,@isnumeric);
p.addParameter('numy',5,@isnumeric);

p.parse(varargin{:});


xlim = p.Results.xlim;
ylim = p.Results.ylim;
sx = p.Results.sx;
sy = p.Results.sy;
safetyFactor = p.Results.safetyFactor;
h0=p.Results.h0;
numx=p.Results.numx;
numy=p.Results.numy;

%% quality investigation
numOnContactGrid = length(find((1-sx)*xlim>mod((xlim+Schnittpunktebene(:,1)),(2*xlim)) |sx*xlim+xlim<mod((xlim+Schnittpunktebene(:,1)),(2*xlim))|(1-sy)*ylim>mod((ylim+Schnittpunktebene(:,2)),(2*ylim)) |sy*ylim+ylim<mod((ylim+Schnittpunktebene(:,2)),(2*ylim))));

%% plot
np=size(Schnittpunktebene,1);
[X,Y] = meshgrid(linspace(-xlim,xlim,numx),linspace(-ylim,ylim,numy));
Z = repmat(-h0,numy,numx);
figure
surface(X,Y,Z)
hold on;
patch([xlim xlim -xlim -xlim sx*xlim sx*xlim xlim],[-ylim ylim ylim sy*ylim sy*ylim -ylim -ylim],zeros(1,7),'FaceColor','k','FaceAlpha',.5);
patch([-xlim -xlim xlim xlim -sx*xlim -sx*xlim -xlim],[ylim -ylim -ylim -sy*ylim -sy*ylim ylim ylim],zeros(1,7),'FaceColor','k','FaceAlpha',.5);
line([orig(:,1)';Schnittpunkt(:,1)';Schnittpunktebene(:,1)'],[orig(:,2)';Schnittpunkt(:,2)';Schnittpunktebene(:,2)'],[orig(:,3)';Schnittpunkt(:,3)';Schnittpunktebene(:,3)'],'Color','k','LineWidth',0.3);
scatter3(Schnittpunktebene(:,1)',Schnittpunktebene(:,2)',zeros(1,np),'*r');
scatter3(sx*safetyFactor*Schnittpunkt(:,1)',sy*safetyFactor*Schnittpunkt(:,2)',zeros(1,np),'.g');
line([Schnittpunktebene(:,1)';sx*safetyFactor*Schnittpunkt(:,1)'], [Schnittpunktebene(:,2)';sy*safetyFactor*Schnittpunkt(:,2)'], zeros(2,np),'Color','r');
title(sprintf('metasurface construction; rays on contact grid: %d/%d',numOnContactGrid,np));
xlabel('x');
ylabel('y');
zlabel('z');
view(0,90)
%axis equal;
hold off


% plot distances
% figure(203);
% dist = reshape(sqrt((Schnittpunktebene(:,1)'-safetyFactor*sx*Schnittpunkte2(:,1)').^2+(Schnittpunktebene(:,2)'-safetyFactor*sy*Schnittpunkte2(:,2)').^2),[],1);
% patch('Faces',tr1.ConnectivityList,'Vertices',tr1.Points(:,1:2),'FaceVertexCData',dist);
% shading flat;
% % pcolor(X,Y,reshape(sqrt((Xhit-sx*Xr).^2+(Yhit-sy*Yr).^2),size(X)));
% % shading interp;
% title(sprintf('distance from design point, average: %2.1f\n# rays on contact grid: %d/%d',mean(sqrt((Schnittpunktebene(:,1)'-safetyFactor*sx*Schnittpunkte2(:,1)').^2+(Schnittpunktebene(:,2)'-safetyFactor*sy*Schnittpunkte2(:,2)').^2)),numOnContactGrid,np));
% colormap jet;
% xlabel('x');
% ylabel('y');
% colorbar;

end

