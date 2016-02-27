% Ray-tracing script for max_refract refractions

%% Parameter
h0=150;
numx=30;
numy=10;
numrayx=50;
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
orig_z=250;
%% angles
alpha=0;
beta=0;
%% construction of the triangulation Object
tr=constructSurface('h0',h0,'numx',numx,'numy',numy,'xlim',xlim,'ylim',ylim,'sx',sx,'sy',sy,'safetyFactor',safetyFactor);
sz=size(tr,1);
%% dir
%constructs the direction of the rays in the beginning
dir=zeros((numrayy-1)*(numrayx-1),3);
dir(:,1)=sind(beta);
dir(:,2)=sind(alpha)*cosd(beta);
dir(:,3)=cosd(alpha)*cosd(beta);
%% orig
%constructs the origin of the rays in the beginning
[x_orig,y_orig]=meshgrid(linspace(xmin+(xlim-xmin)/(4*numrayx-4),xlim-(xlim-xmin)/(4*numrayx-4),numrayx-1),linspace(ymin+(ylim-ymin)/(4*numrayy-4),ylim-(ylim-ymin)/(4*numrayy-4),numrayy-1));
z_orig=zeros((numrayy-1)*(numrayx-1),1);
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
triangular_Point1=Points(Connect(:,1),:);
triangular_Point2=Points(Connect(:,2),:);
triangular_Point3=Points(Connect(:,3),:);
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
Schnittpunkte_aktuell = squeeze(Schnittpunkte(k-1,:,:));
vtrans_aktuell = squeeze(vtrans(k-1,:,:));
% [vtrans_aktuell,~] = cloneToEqualSize(vtrans_aktuell,triangular_Point1);
% [triangular_Point1,~] = cloneToEqualSize(triangular_Point1,Schnittpunkte_aktuell);
% [triangular_Point2,~] = cloneToEqualSize(triangular_Point2,Schnittpunkte_aktuell);
% [triangular_Point3,Schnittpunkte_aktuell] = cloneToEqualSize(triangular_Point3,Schnittpunkte_aktuell);

[triInd,rayInd] = meshgrid(1:size(triangular_Point1,1),1:size(Schnittpunkte_aktuell,1));
tic;
intersect = TriangleRayIntersection(Schnittpunkte_aktuell(rayInd(:),:),vtrans_aktuell(rayInd(:),:),triangular_Point1(triInd(:),:),triangular_Point2(triInd(:),:),triangular_Point3(triInd(:),:));
intersect = reshape(intersect,size(triInd));
toc
tic;
[Schnittpunkte(k,:,:),Normale]=calculateSPofRayswithSurface(tr,triangular_Point1,triangular_Point2,triangular_Point3,reshape(Schnittpunkte(k-1,:,:),size(Schnittpunkte,2),size(Schnittpunkte,3)),reshape(vtrans(k-1,:,:),size(vtrans,2),size(vtrans,3)));
toc
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


%% plotting

plotSurface(tr,Schnittpunkte,Schnittpunktebene,k,'sx',sx,'sy',sy,'safetyFactor',safetyFactor);

