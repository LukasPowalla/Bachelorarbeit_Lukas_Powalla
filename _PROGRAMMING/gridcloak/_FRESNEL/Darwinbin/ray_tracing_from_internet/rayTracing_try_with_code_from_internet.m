% Ray-tracing script

%% Parameter
h0=30;
numx=20;
numy=5;
numrayx=50;
numrayy=5;
xlim=500;
ylim=50;
%% orig
tr=constructSurface('h0',h0,'numx',numx,'numy',numy,'xlim',xlim,'ylim',ylim);
sz=size(tr,1);
[x_orig,y_orig]=meshgrid(linspace(0,xlim,numrayx),linspace(0,ylim,numrayy));
z_orig=zeros(numrayy*numrayx,1)-500;
orig=[x_orig(:),y_orig(:),z_orig];
%% dir
dir=zeros(numrayy*numrayx,1);
dir(:,3)=1;
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
orig=repmat(orig,size(triangular_Point1,1),1);
dir=repmat(dir,size(sz,1),1);
TriangleRayIntersection(orig, dir, triangular_Point1, triangular_Point2,triangular_Point3);