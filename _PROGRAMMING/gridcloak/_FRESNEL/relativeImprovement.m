function relimp = relativeImprovement(tr,alpha,beta,varargin)
%% plot optimized surface and raytracing results

%% input parsing
p = inputParser;

p.addParameter('n1',1,@isnumeric);              % starting ref. index
p.addParameter('n2',1.5,@isnumeric);            % end ref. index
p.addParameter('xlim',500,@isnumeric);          % x-size of unit cell
p.addParameter('ylim',50,@isnumeric);           % y-size of unit cell
p.addParameter('sx',0.9,@isnumeric);            % x-scaling factor
p.addParameter('sy',0.9,@isnumeric);            % y-scaling factor
p.addParameter('safetyFactor',0.9,@isnumeric);  % additional safety factor
p.addParameter('h0',100,@isnumeric);  
p.addParameter('orig_z',200,@isnumeric);  
p.addParameter('numrayx',50,@isnumeric);  
p.addParameter('numrayy',10,@isnumeric);  
p.addParameter('xmin',-500,@isnumeric);  
p.addParameter('ymin',-50,@isnumeric);  

p.parse(varargin{:});

n1 = p.Results.n1;
n2 = p.Results.n2;
xlim = p.Results.xlim;
ylim = p.Results.ylim;
sx = p.Results.sx;
sy = p.Results.sy;
safetyFactor=p.Results.safetyFactor;
h0=p.Results.h0;
numrayx=p.Results.numrayx;
numrayy=p.Results.numrayy;
xmin=p.Results.xmin;
ymin=p.Results.ymin;
orig_z=p.Results.orig_z;
%% size of triangulation Object
sz=size(tr,1);
%% dir
%constructs the direction of the rays in the beginning
dir=zeros((numrayy)*(numrayx),3);
dir(:,1)=sind(beta);
dir(:,2)=sind(alpha)*cosd(beta);
dir(:,3)=cosd(alpha)*cosd(beta);
%% orig
%constructs the origin of the rays in the beginning
[x_orig,y_orig]=meshgrid(linspace(xmin+(xlim-xmin)/(4*numrayx),xlim-(xlim-xmin)/(4*numrayx),numrayx),linspace(ymin+(ylim-ymin)/(4*numrayy),ylim-(ylim-ymin)/(4*numrayy),numrayy));
z_orig=zeros((numrayy)*(numrayx),1);
orig=[x_orig(:),y_orig(:),z_orig];
orig=(-(orig_z-h0)/(cosd(alpha)*cosd(beta)))*dir+orig;
orig(:,3)=orig(:,3)-h0;
%% vert 
% constructs an object conststing out of
% triangular_Point1,triangular_Point2,triangular_Point3 describing the
% edges of every triangular needed for the calculateSPofRayswithSurface
% function
Points=tr.Points;
Connect=tr.ConnectivityList;
triangular_Point1=zeros(sz,3);
triangular_Point2=zeros(sz,3);
triangular_Point3=zeros(sz,3);
for i=1:sz
    triangular_Point1(i,:)=Points(Connect(i,1),:);
    triangular_Point2(i,:)=Points(Connect(i,2),:);
    triangular_Point3(i,:)=Points(Connect(i,3),:);
end
%% calculates refraction for maximum max_refract refractions
%defining all the arrays needed to know 
%Schnittpunkte are the Intersection points with the Triangulation Object
%vtrans is the array containing all the velocityvectors
k=1;
max_refract=20;
Schnittpunkte=nan(max_refract,size(orig,1),3);
vtrans=nan(max_refract,size(orig,1),3);
Schnittpunkte(1,:,:)=orig;
vtrans(1,:,:)=dir;
% calculate recraction until all of the rays don't hit the triangulation
% Object any more or until max_refract is reached
while ~isnan(min((Schnittpunkte(k,:,3))))&&k<max_refract
k=k+1;
[Schnittpunkte(k,:,:),Normale]=calculateSPofRayswithSurface(tr,triangular_Point1,triangular_Point2,triangular_Point3,reshape(Schnittpunkte(k-1,:,:),size(Schnittpunkte,2),size(Schnittpunkte,3)),reshape(vtrans(k-1,:,:),size(vtrans,2),size(vtrans,3)));
vtrans(k,:,:)=Snell(reshape(vtrans(k-1,:,:),size(vtrans,2),size(vtrans,3)),Normale,n1,n2);
end
%For plotting it is useful to copy the values to the next array in case
%they don't hit the triangulation object any more. This is useful because
%than you have in Schnittpunkte(k-1) all the endpoints of the refraction
%and in vtrans(k-1) all the endvelocities
for m=2:(k-1)
    for l=1:size(Schnittpunkte,2)
        for j=1:3
            if isnan(Schnittpunkte(m,l,j))
                Schnittpunkte(m,l,j)=Schnittpunkte(m-1,l,j);
                vtrans(m,l,j)=vtrans(m-1,l,j);
            end
        end
    end
end
%now, you calculate the Intersection with the plane
Schnittpunktebene=zeros(length(orig),3);
Schnittpunktebene(:,1)=Schnittpunkte((k-1),:,1)'-vtrans((k-1),:,1)'.*(vtrans((k-1),:,3)'.\Schnittpunkte((k-1),:,3)');
Schnittpunktebene(:,2)=Schnittpunkte((k-1),:,2)'-vtrans((k-1),:,2)'.*(vtrans((k-1),:,3)'.\Schnittpunkte((k-1),:,3)');
Schnittpunktebene(:,3)=Schnittpunkte((k-1),:,3)'-vtrans((k-1),:,3)'.*(vtrans((k-1),:,3)'.\Schnittpunkte((k-1),:,3)');

numOnContactGrid = length(find((1-sx)*xlim>mod((xlim+Schnittpunktebene(:,1)),(2*xlim)) |sx*xlim+xlim<mod((xlim+Schnittpunktebene(:,1)),(2*xlim))|(1-sy)*ylim>mod((ylim+Schnittpunktebene(:,2)),(2*ylim)) |sy*ylim+ylim<mod((ylim+Schnittpunktebene(:,2)),(2*ylim))));
np=size(Schnittpunktebene,1);
%% calculation of the relative Improvement
f=1-sx*sy;
N0=np;
N=np-numOnContactGrid;
relimp=N/(N0*(1-f))-1;
end