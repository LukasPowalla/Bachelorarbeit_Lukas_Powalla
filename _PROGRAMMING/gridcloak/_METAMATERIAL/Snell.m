function [ v_trans ] = Snell( dphix,dphiy,v,n,alpha,beta,varargin)
%generalised snell's law
%  normal vectors in direction of v 
% v and n need to be normalised
% v and n are n x 3
%% input parsing
p = inputParser;

p.addParameter('n1',1,@isnumeric);              % starting ref. index
p.addParameter('n2',1.5,@isnumeric);            % end ref. index
p.addParameter('k0',1,@isnumeric);

p.parse(varargin{:});

n1 = p.Results.n1;
n2 = p.Results.n2;
k0=p.Results.k0;
%% preparing angles
teta_i=cosd(v*n');
gamma=atand(sind(alpha)./tand(beta));
if beta==0
   gamma=90; 
end
if alpha==0
   gamma=0; 
end
%% preparing the new dphi derivatives
dphiv=-sind(gamma).*dphix+cosd(gamma).*dphiy;
dphiu=cosd(gamma).*dphix+sind(gamma).*dphiy;
%% calculation of the angulars of the refracted rays in the coordinate system u/v 
% u=cos(gamma)*x+sin(gamma)*y
% v=-sin(gamma)*x+cos(gamma)*y
theta_t=asind((dphiu./k0+n1.*sind(teta_i))./n2);
phi_t=asind(dphiv./(n2.*k0.*cosd(theta_t)));
%% transformation back to x,y and z coordinates
v_trans(:,1)=tand(theta_t(:));
v_trans(:,2)=sind(phi_t(:));
v_trans(:,3)=cosd(phi_t(:));
for i=1:size(dphix,1)*size(dphix,2)
    v_trans(i,:)=v_trans(i,:)/norm(v_trans(i,:));
end

