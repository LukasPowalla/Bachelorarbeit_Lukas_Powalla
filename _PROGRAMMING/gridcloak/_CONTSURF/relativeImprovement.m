function relimp = relativeImprovement(X,Y,Z,alpha,beta,varargin)
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

p.parse(varargin{:});

n1 = p.Results.n1;
n2 = p.Results.n2;
xlim = p.Results.xlim;
ylim = p.Results.ylim;
sx = p.Results.sx;
sy = p.Results.sy;

%% convert meshgrid to triangulated surface
[P,N,tr] = makeTri(X,Y,Z);
v0 = [sind(beta) sind(alpha).*cosd(beta) cosd(alpha).*cosd(beta)];
 Winkel=acosd(N*v0');
 ExcludeIndex=(Winkel>=90);
% Alle Zeilennummer mit Zeilennummer in ExcludeIndex müssen aus Objekt
% [P,N,tr] entfernt werden
P = P(~ExcludeIndex,:);
N = N(~ExcludeIndex,:);
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
%% calculation of the relative Improvement
f=1-sx*sy;
N0=np;
N=np-numOnContactGrid;
relimp=N/(N0*(1-f))-1;
end