function [ v_trans ] = Snell( dphix,dphiy,dxz,alpha,beta,varargin)
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
theta_i= -acosd((cosd(alpha).*cosd(beta)-dxz.*sind(beta))./sqrt(dxz.^2+1));
dphiu=-sqrt((dxz.^2+1)).*cosd(beta).*sqrt(1+((dxz.^2)-cosd(alpha).^2).*cosd(beta).^2+2.*cosd(alpha).*cosd(beta).*sind(beta).*dxz).*sind(alpha)./(((dxz.^2).*cosd(alpha).^2+sind(alpha).^2.*(dxz.^2+1)).*cosd(beta).^2+2.*cosd(alpha).*cosd(beta).*sind(beta).*dxz+sind(beta).^2).*dphiy;
dphiv=sqrt(1+(dxz.^2-cosd(alpha).^2).*cosd(beta).^2+2.*cosd(alpha).*cosd(beta).*sind(beta).*dxz).*(cosd(alpha).*cosd(beta).*dxz+sind(beta))./((dxz.^2.*cosd(alpha).^2+sind(alpha).^2.*(dxz.^2+1)).*cosd(beta).^2+2.*cosd(alpha).*cosd(beta).*sind(beta).*dxz+sind(beta).^2).*dphiy;
% dphiw=0;

theta_t=asind((dphiu./k0+n1.*sind(theta_i))./n2);
phi_t=asind(dphiv./(n2.*k0.*cosd(theta_t)));
%% transformation back to x,y and z coordinates

v_xyz_snell_x=-(sind(alpha).*cosd(beta).*cosd(theta_t).*sind(phi_t).*sqrt((dxz.^2+1))+dxz.*cosd(phi_t).*cosd(theta_t).*sqrt(-cosd(alpha).^2.*cosd(beta).^2+2.*cosd(alpha).*cosd(beta).*sind(beta).*dxz+cosd(beta).^2.*(dxz.^2)+1)+cosd(alpha).*cosd(beta).*sind(theta_t).*dxz+sind(beta).*sind(theta_t)).*((dxz.^2+1).^(-1./2)).*(-cosd(alpha).^2.*cosd(beta).^2+2.*cosd(alpha).*cosd(beta).*sind(beta).*dxz+cosd(beta).^2.*(dxz.^2)+1).^(-1./2);
v_xyz_snell_y=-(-cosd(alpha).*cosd(beta).*cosd(theta_t).*sind(phi_t).*dxz+sqrt(dxz.^2+1).*sind(alpha).*cosd(beta).*sind(theta_t)-sind(beta).*cosd(theta_t).*sind(phi_t)).*(-cosd(alpha).^2.*cosd(beta).^2+2.*cosd(alpha).*cosd(beta).*sind(beta).*dxz+cosd(beta).^2.*dxz.^2+1).^(-1./2);
v_xyz_snell_z=-(dxz.*sind(alpha).*cosd(beta).*cosd(theta_t).*sind(phi_t).*sqrt((dxz.^2+1))+cosd(alpha).*cosd(beta).*sind(theta_t).*(dxz.^2)-cosd(phi_t).*cosd(theta_t).*sqrt(-cosd(alpha).^2.*cosd(beta).^2+2.*cosd(alpha).*cosd(beta).*sind(beta).*dxz+cosd(beta).^2.*(dxz.^2)+1)+sind(beta).*sind(theta_t).*dxz).*((dxz.^2+1).^(-1./2)).*(-cosd(alpha).^2.*cosd(beta).^2+2.*cosd(alpha).*cosd(beta).*sind(beta).*dxz+cosd(beta).^2.*(dxz.^2)+1).^(-1./2);
v_trans=[v_xyz_snell_x(:) v_xyz_snell_y(:) v_xyz_snell_z(:) ];

