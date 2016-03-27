function [ dphix,dphiy,x_m,y_m ] = constructSurface( varargin )
%constructs Surface, especially the matrices dphix and dphiy
%% input parsing
p = inputParser;

p.addParameter('n1',1,@isnumeric);              % starting ref. index
p.addParameter('n2',1.5,@isnumeric);            % end ref. index
p.addParameter('xlim',500,@isnumeric);          % x-size of unit cell
p.addParameter('ylim',50,@isnumeric);           % y-size of unit cell
p.addParameter('numx',50,@isnumeric);
p.addParameter('numy',10,@isnumeric);
p.addParameter('sx',0.9,@isnumeric);            % x-scaling factor
p.addParameter('sy',0.9,@isnumeric);            % y-scaling factor
p.addParameter('safetyFactor',0.9,@isnumeric);  % additional safety factor
p.addParameter('h0',20,@isnumeric);
p.addParameter('k0',1,@isnumeric);


p.parse(varargin{:});


n1 = p.Results.n1;
n2 = p.Results.n2;
xlim = p.Results.xlim;
ylim = p.Results.ylim;
sx = p.Results.sx;
sy = p.Results.sy;
safetyFactor = p.Results.safetyFactor;
numx=p.Results.numx;
numy=p.Results.numy;
h0=p.Results.h0;
k0=p.Results.k0;
%% angular preparation
alpha=0;
beta=0;
gamma=atand(sind(alpha)/tand(beta));
if beta==0
   gamma=90;
end
if alpha==0
   gamma=0; 
end
theta_i=acosd(cosd(alpha)*cosd(beta));
%theta_i=atand((-cosd(gamma)*sind(beta)-sind(gamma)*sind(alpha)*cosd(beta))/(-cosd(alpha)*cosd(beta)));
%%
xmin=-xlim;
ymin=-ylim;
sy=sy*safetyFactor;
sx=sx*safetyFactor;
%% x_m and y_m are the midpoints of the pixels, we have numx*numy pixels
[x_m,y_m]=meshgrid(linspace(xmin+(xlim-xmin)/(2*numx),xlim-(xlim-xmin)/(2*numx),numx),linspace(ymin+(ylim-ymin)/(2*numy),ylim-(ylim-ymin)/(2*numy),numy));
phi_t=atand((-(sx-1).*x_m.*sind(gamma)+cosd(gamma).*(sy-1).*y_m)./h0);
arg1=-(sx-1).*x_m.*cosd(gamma)-(sy-1).*y_m.*sind(gamma);
arg2=sqrt((-(sx-1).*x_m.*sind(gamma)+cosd(gamma).*(sy-1).*y_m).^2+h0.^2);
theta_t=atand(arg1./arg2);
dphiu=(k0.*(n2.*sind(theta_t)-n1.*sind(theta_i)));
dphiv=n2.*k0.*cosd(theta_t).*sind(phi_t);
dphix=-cosd(gamma).*dphiu-sind(gamma).*dphiv;
dphiy=-sind(gamma).*dphiu+cosd(gamma).*dphiv;
% u_m=-cosd(gamma).*x_m-sind(gamma).*y_m;
% v_m=-sind(gamma).*x_m+cosd(gamma).*y_m;
end

