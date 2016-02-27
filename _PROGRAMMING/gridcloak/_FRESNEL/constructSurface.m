function [tr1] = constructSurface(varargin)
% constructing fresnel optic
%% input parsing
p = inputParser;

p.addParameter('n1',1,@isnumeric);              % starting ref. index
p.addParameter('n2',1.5,@isnumeric);            % end ref. index
p.addParameter('xlim',500,@isnumeric);          % x-size of unit cell
p.addParameter('ylim',50,@isnumeric);           % y-size of unit cell
p.addParameter('sx',0.9,@isnumeric);            % x-scaling factor
p.addParameter('sy',0.9,@isnumeric);            % y-scaling factor
p.addParameter('safetyFactor',0.9,@isnumeric);  % additional safety factor
p.addParameter('h0',150,@isnumeric);
p.addParameter('numx',31,@isnumeric);
p.addParameter('numy',4,@isnumeric);

p.parse(varargin{:});

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


%% Parameters
sx=sx*safetyFactor;
sy=sy*safetyFactor;
xmin=0;
ymin=0;
%% meshgrids
[x,y]=meshgrid(linspace(xmin,xlim,numx),linspace(ymin,ylim,numy));
[x_m,y_m]=meshgrid(linspace(xmin+(xlim-xmin)/(2*numx-2),xlim-(xlim-xmin)/(2*numx-2),numx-1),linspace(ymin+(ylim-ymin)/(2*numy-2),ylim-(ylim-ymin)/(2*numy-2),numy-1));
[nx,ny,nz]=calculateNormalvector(x_m,y_m,'h0',h0,'n1',n1,'n2',n2,'sx',sx,'sy',sy);
%% calculating Points for triangulation
% we are calculating all the points of the surface at x,y using the
% calculated normal vectors
Points=zeros(4*(numx-1)*(numy-1)+2*(numy-1),3);
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
% dealing with the edges
m=(numx-1)*(numy-1)*4+1;
ye=linspace(ymin,ylim,numy);
for i=1:(numy-1)
    Points(m,:)=[xlim ye(i) -h0];
    Points(m+1,:)=[xlim ye(i+1) -h0];
    m=m+2;
end


%% calculating connectivity-list for triangulation
ConnectivityList=[];
% construct triangles of "pixels"
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
% constructing triangles of vertical surfaces in x- direction
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
% constructing triangles of vertical surfaces in the other direction
m=(2*(numy-1)*(numx-2))+2*(numy-1)*(numx-1)+1;
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
% dealing with the edges
for i=0:(numy-2)
    ConnectivityList(m,1)=(numx-1)*(numy-1)*4+1+2*i;
    ConnectivityList(m,2)=(numx-1)*(numy-1)*4+2+2*i;
    ConnectivityList(m,3)=4*(numx-1)-1+4*i*(numx-1);
    ConnectivityList(m+1,1)=4*(numx-1)-2+4*i*(numx-1);
    ConnectivityList(m+1,2)=(numx-1)*(numy-1)*4+1+2*i;
    ConnectivityList(m+1,3)=4*(numx-1)-1+4*i*(numx-1);
    m=m+2;
end

numberofpoints=size(Points(:,1),1);
Points1=[-Points(:,1) Points(:,2) Points(:,3) ];
Points2=[Points(:,1) -Points(:,2) Points(:,3) ];
Points3=[-Points(:,1) -Points(:,2) Points(:,3) ];
%% triangulation object
Points_gesamt=[Points;Points1;Points2;Points3];
ConnectivityList_gesamt=[ConnectivityList;(ConnectivityList+numberofpoints);(ConnectivityList+2*numberofpoints);(ConnectivityList+3*numberofpoints)];
tr1=triangulation(ConnectivityList_gesamt,Points_gesamt);
% plot surface
figure(132)
trisurf(tr1,'FaceColor','interp');
hold on;
patch([xlim xlim -xlim -xlim sx*xlim sx*xlim xlim],[-ylim ylim ylim sy*ylim sy*ylim -ylim -ylim],zeros(1,7),'FaceColor','k','FaceAlpha',.5);
patch([-xlim -xlim xlim xlim -sx*xlim -sx*xlim -xlim],[ylim -ylim -ylim -sy*ylim -sy*ylim ylim ylim],zeros(1,7),'FaceColor','k','FaceAlpha',.5);
title('fresnel construction');
xlabel('x');
ylabel('y');
zlabel('z');
%axis equal;
hold off

end

