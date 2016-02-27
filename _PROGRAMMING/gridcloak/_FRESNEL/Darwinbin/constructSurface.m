% constructing fresnel optic
%% Parameters
sx=0.9;
sy=0.9;
numx=20;
numy=5;
xmin=0;
xmax=500;
ymin=0;
ymax=50;
h0=25;

[x,y]=meshgrid(linspace(xmin,xmax,numx),linspace(ymin,ymax,numy));
[x_m,y_m]=meshgrid(linspace(xmin+(xmax-xmin)/(2*numx-2),xmax-(xmax-xmin)/(2*numx-2),numx-1),linspace(ymin+(ymax-ymin)/(2*numy-2),ymax-(ymax-ymin)/(2*numy-2),numy-1));
[nx,ny,nz]=calculateNormalvector(x_m,y_m,'h0',h0);
%% calculating Points for triangulation
Points=zeros(4*(numx-1)*(numy-1),3);
m=1;
for i=1:(numy-1)
    for k=1:(numx-1)
        z1=(nx(i,k)*(x(1,k)-x_m(i,k))+ny(i,k)*(y(i,1)-y_m(i,k)))/nz(i,k)-h0;
        z2=(nx(i,k)*(x(1,k+1)-x_m(i,k))+ny(i,k)*(y(i,1)-y_m(i,k)))/nz(i,k)-h0;
        z3=(nx(i,k)*(x(1,k+1)-x_m(i,k))+ny(i,k)*(y(i+1,1)-y_m(i,k)))/nz(i,k)-h0;
        z4=(nx(i,k)*(x(1,k)-x_m(i,k))+ny(i,k)*(y(i+1,1)-y_m(i,k)))/nz(i,k)-h0;
        Points(m,:)=[x(1,k) y(i,1) z1];
        Points(m+1,:)=[x(1,k+1) y(i,1) z2];
        Points(m+2,:)=[x(1,k+1) y(i+1,1) z3];
        Points(m+3,:)=[x(1,k) y(i+1,1) z4];
        m=m+4;
    end
end
%% calculating connectivity-list for triangulation
ConnectivityList=[];%
m=1;
for i=0:((numy-1)*(numx-1)-1)
    ConnectivityList(m,1)=1+4*i;
    ConnectivityList(m,2)=2+4*i;
    ConnectivityList(m,3)=4+4*i;
    ConnectivityList(m+1,1)=2+4*i;
    ConnectivityList(m+1,2)=3+4*i;
    ConnectivityList(m+1,3)=4+4*i;
    m=m+2;
end
m=(2*(numy-1)*(numx-1)+1);
for k=0:(numy-2)
    for i=0:(numx-3)
        ConnectivityList(m,1)=2+4*i+4*(numx-1)*k;
        ConnectivityList(m,2)=3+4*i+4*(numx-1)*k;
        ConnectivityList(m,3)=5+4*i+4*(numx-1)*k;
        ConnectivityList(m+1,1)=3+4*i+4*(numx-1)*k;
        ConnectivityList(m+1,2)=5+4*i+4*(numx-1)*k;
        ConnectivityList(m+1,3)=8+4*i+4*(numx-1)*k;
        m=m+2;
    end
end
m=(2*(numy-1)*(numx-1)+1)+2*(numy-1)*(numx-2);
for k=0:(numy-3)
    for i=0:(numx-2)
        ConnectivityList(m,1)=3+4*i+k*4*(numx-1);
        ConnectivityList(m,2)=4+4*i+k*4*(numx-1);
        ConnectivityList(m,3)=1+4*(numx-1)+4*i+k*4*(numx-1);
        ConnectivityList(m+1,1)=2+4*(numx-1)+4*i+k*4*(numx-1);
        ConnectivityList(m+1,2)=3+4*i+k*4*(numx-1);
        ConnectivityList(m+1,3)=1+4*(numx-1)+4*i+k*4*(numx-1);
        m=m+2;
    end
end

numberofpoints=size(Points(:,1),1);
Points1=[-Points(:,1) Points(:,2) Points(:,3) ];
Points2=[Points(:,1) -Points(:,2) Points(:,3) ];
Points3=[-Points(:,1) -Points(:,2) Points(:,3) ];
%% triangulation object
Points_gesamt=[Points;Points1;Points2;Points3];
ConnectivityList_gesamt=[ConnectivityList;(ConnectivityList+numberofpoints);(ConnectivityList+2*numberofpoints);(ConnectivityList+3*numberofpoints)];
tr1=triangulation(ConnectivityList_gesamt,Points_gesamt);
%% plot surface
figure(132)
trisurf(tr1,'FaceColor','interp');
hold on;
patch([xmax xmax -xmax -xmax sx*xmax sx*xmax xmax],[-ymax ymax ymax sy*ymax sy*ymax -ymax -ymax],zeros(1,7),'FaceColor','k','FaceAlpha',.5);
patch([-xmax -xmax xmax xmax -sx*xmax -sx*xmax -xmax],[ymax -ymax -ymax -sy*ymax -sy*ymax ymax ymax],zeros(1,7),'FaceColor','k','FaceAlpha',.5);
title('fresnel construction');
xlabel('x');
ylabel('y');
zlabel('z');
%axis equal;
hold off


