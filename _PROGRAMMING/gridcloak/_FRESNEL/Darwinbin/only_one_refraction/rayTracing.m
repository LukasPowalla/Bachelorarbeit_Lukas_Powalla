% Ray-tracing script

%% Parameter
h0=80;
numx=20;
numy=5;
numrayx=20;
numrayy=5;
xlim=500;
ylim=50;
xmin=-500;
ymin=-50;
n1=1;
n2=1.5;
safetyFactor=0.9;
sx=0.9;
sy=0.9;
orig_z=200;
%% angles
alpha=0;
beta=0;
%% dir
dir=zeros((numrayy-1)*(numrayx-1),3);
dir(:,1)=sind(beta);
dir(:,2)=sind(alpha)*cosd(beta);
dir(:,3)=cosd(alpha)*cosd(beta);
%% orig
tr=constructSurface('h0',h0,'numx',numx,'numy',numy,'xlim',xlim,'ylim',ylim,'sx',sx,'sy',sy);
sz=size(tr,1);
[x_orig,y_orig]=meshgrid(linspace(xmin+(xlim-xmin)/(4*numrayx-4),xlim-(xlim-xmin)/(4*numrayx-4),numrayx-1),linspace(ymin+(ylim-ymin)/(4*numrayy-4),ylim-(ylim-ymin)/(4*numrayy-4),numrayy-1));
z_orig=zeros((numrayy-1)*(numrayx-1),1);
orig=[x_orig(:),y_orig(:),z_orig];
orig=(-orig_z/(cosd(alpha)*cosd(beta)))*dir+orig;
%% vert i
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

[Schnittpunkte1,Normale1]=calculateSPofRayswithSurface(tr,triangular_Point1,triangular_Point2,triangular_Point3,orig,dir);
vtrans1=Snell(dir,Normale1,n1,n2);
%Schnittpunktebene=zeros(length(orig),3);
% Schnittpunktebene(:,1)=Schnittpunkte1(:,1)-vtrans1(:,1).*(vtrans1(:,3).\Schnittpunkte1(:,3));
% Schnittpunktebene(:,2)=Schnittpunkte1(:,2)-vtrans1(:,2).*(vtrans1(:,3).\Schnittpunkte1(:,3));
% Schnittpunktebene(:,3)=Schnittpunkte1(:,3)-vtrans1(:,3).*(vtrans1(:,3).\Schnittpunkte1(:,3));

[Schnittpunkte2,Normale2]=calculateSPofRayswithSurface(tr,triangular_Point1,triangular_Point2,triangular_Point3,Schnittpunkte1,vtrans1);
vtrans2=Snell(vtrans1,Normale2,n1,n2);

for i=1:size(dir,1)
    if isnan(Schnittpunkte2(i,:))
       Schnittpunkte2(i,:)=Schnittpunkte1(i,:); 
    end
    if isnan(vtrans2(i,:))
       vtrans2(i,:)=vtrans1(i,:); 
    end
end

Schnittpunktebene=zeros(length(orig),3);
Schnittpunktebene(:,1)=Schnittpunkte2(:,1)-vtrans2(:,1).*(Schnittpunkte2(:,3)./vtrans2(:,3));
Schnittpunktebene(:,2)=Schnittpunkte2(:,2)-vtrans2(:,2).*(Schnittpunkte2(:,3)./vtrans2(:,3));
Schnittpunktebene(:,3)=Schnittpunkte2(:,3)-vtrans2(:,3).*(Schnittpunkte2(:,3)./vtrans2(:,3));
%% plotting
%Schnittpunkte2=Schnittpunkte1;
plotSurface(tr,orig,Schnittpunkte1,Schnittpunkte2,Schnittpunktebene,'sx',sx,'sy',sy);

