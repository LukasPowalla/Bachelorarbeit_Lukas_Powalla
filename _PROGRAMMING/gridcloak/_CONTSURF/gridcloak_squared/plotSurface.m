function fitness = plotSurface(X,Y,Z,varargin)
%% plot optimized surface and raytracing results

%% input parsing
p = inputParser;

p.addOptional('alpha',0,@isnumeric);
p.addOptional('beta',0,@isnumeric);
p.addParameter('n1',1,@isnumeric);              % starting ref. index
p.addParameter('n2',1.5,@isnumeric);            % end ref. index
p.addParameter('xlim',50,@isnumeric);          % x-size of unit cell
p.addParameter('ylim',50,@isnumeric);           % y-size of unit cell
p.addParameter('sx',0.9,@isnumeric);            % x-scaling factor
p.addParameter('sy',0.9,@isnumeric);            % y-scaling factor
p.addParameter('safetyFactor',0.9,@isnumeric);  % additional safety factor

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

%% convert meshgrid to triangulated surface
[P,N,tr] = makeTri(X,Y,Z);
v0 = [sind(beta) sind(alpha)*cosd(beta) cosd(alpha)*cosd(beta)];
 Winkel=acosd(N*v0');
 ExcludeIndex=(Winkel>=90);
% Alle Zeilennummer mit Zeilennummer in ExcludeIndex müssen aus Objekt
% [P,N,tr] entfernt werden
tr0 = tr;
P = P(~ExcludeIndex,:);
N = N(~ExcludeIndex,:);
trconnlist = tr.ConnectivityList(~ExcludeIndex,:);
% trpoints = tr.Points(trconnlist(ismember(1:length(tr.Points),trconnlist(:))),:);

warnstate = warning;
warning('off','MATLAB:triangulation:PtsNotInTriWarnId');
tr = triangulation(trconnlist,tr.Points);
warning(warnstate);

np = size(P,1);


%% raytracing
v0 = repmat(v0,np,1);
Xhit = nan*ones(np,1);
Yhit = nan*ones(np,1);

v1=Snell(v0,N,n1,n2);

z0 = min(min(Z))*1.1;
A=repmat([-1 0; 0 -1; 0 0],1,1,np);
A(:,3,:)=reshape(v1',3,1,[]);

Xr = P(:,1);
Yr = P(:,2);
Zr = P(:,3);
B = [-Xr, -Yr, -Zr];
B = reshape(B',3,1,[]);

for m = 1:size(B,3)
    res = A(:,:,m)\B(:,:,m);
    Xhit(m) = res(1);
    Yhit(m) = res(2);
end

% count rays hitting contact finger
%numOnContactGrid1 = length(find(Xhit>sx*xlim | Yhit>sy*ylim | Xhit<-sx*xlim | Yhit<-sy*ylim));
numOnContactGrid = length(find((1-sx)*xlim>mod((xlim+Xhit),(2*xlim)) |sx*xlim+xlim<mod((xlim+Xhit),(2*xlim))|(1-sy)*ylim>mod((ylim+Yhit),(2*ylim)) |sy*ylim+ylim<mod((ylim+Yhit),(2*ylim))));
% test=mod((xlim+Xhit),(2*xlim));
% test2=sx*xlim;
%% calculate starting points in dependance of angle
A=repmat([1 0; 0 1; 0 0],1,1,np);
A(:,3,:)=reshape(v0',3,1,[]);

Xr = P(:,1);
Yr = P(:,2);
Zr = P(:,3);
B = [Xr, Yr, Zr-z0];
B = reshape(B',3,1,[]);

Px = nan*ones(np,1);
Py = nan*ones(np,1);

for m = 1:size(B,3)
    res = A(:,:,m)\B(:,:,m);
    Px(m) = res(1);
    Py(m) = res(2);
end


%% plotting
if nargout==0
% plot surface
figure(201);
trisurf(tr0);
title('surface height');
xlabel('x');
ylabel('y');
zlabel('z');
colorbar;
% caxis([min(Zr),max(Zr)]);
axis equal;
hold off;

% plot rays
figure(202);
trisurf(tr0);
hold on;
patch([xlim xlim -xlim -xlim sx*xlim sx*xlim xlim],[-ylim ylim ylim sy*ylim sy*ylim -ylim -ylim],zeros(1,7),'FaceColor','k','FaceAlpha',.5);
patch([-xlim -xlim xlim xlim -sx*xlim -sx*xlim -xlim],[ylim -ylim -ylim -sy*ylim -sy*ylim ylim ylim],zeros(1,7),'FaceColor','k','FaceAlpha',.5);
line([Px'; Xr'; Xhit'],[Py'; Yr'; Yhit'],[repmat(z0,1,np); Zr'; zeros(1,np)],'Color','k','LineWidth',0.3);
scatter3(Xhit,Yhit,zeros(1,np),'*r');
scatter3(safetyFactor*sx*Xr,safetyFactor*sy*Yr,zeros(1,np),'.g');
line([Xhit';safetyFactor*sx*Xr'], [Yhit';safetyFactor*sy*Yr'], zeros(2,np),'Color','r');
title('ray tracing');
xlabel('x');
ylabel('y');
zlabel('z');
colorbar;
view(2);
% caxis([min(Zr),max(Zr)]);
hold off;


% plot distances
figure(203);
dist = reshape(sqrt((Xhit-safetyFactor*sx*Xr).^2+(Yhit-safetyFactor*sy*Yr).^2),[],1);
patch('Faces',tr.ConnectivityList,'Vertices',tr.Points(:,1:2),'FaceVertexCData',dist);
shading flat;
% pcolor(X,Y,reshape(sqrt((Xhit-sx*Xr).^2+(Yhit-sy*Yr).^2),size(X)));
% shading interp;
title(sprintf('distance from design point, average: %2.1f\n# rays on contact grid: %d/%d',mean(sqrt((Xhit-safetyFactor*sx*Xr).^2+(Yhit-safetyFactor*sy*Yr).^2)),numOnContactGrid,np));
colormap jet;
xlabel('x');
ylabel('y');
colorbar;

else

fitness = mean(sqrt((Xhit-safetyFactor*sx*Xr).^2+(Yhit-safetyFactor*sy*Yr).^2)) * (numOnContactGrid+1);

end
clear xlim ylim
end