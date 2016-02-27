function [TavgCloak,TavgBare,TavgRef] = FresnelDif(P,Pref,n1,n2,reflections)
%% returns the angle-averaged transmittance for diffusive irradiance
% OUTPUTS
%  TavgCloak: averaged transmittance with cloak
%  TavgBare: averaged transmittance without any structure
%  TavgRef: averaged transmittance with reference structure
% INPUTS
%  P: interpolant giving P(alpha,beta) with cloak
%  Pref: interpolant giving P(alpha,beta) with reference structure
%  s: settings container

%% parameters
% angle discretization
da = 0.5;

%% input parsing
maxs = max(P.Points);
mins = min(P.Points);
alphas = mins(1):da:maxs(1);
betas = mins(2):da:maxs(2);

if isempty(Pref)
    referenced = false;
else
    referenced = true;
end

%% coordinates
% modified spherical coordinates (as used in rayTracer):
%  y-axis is direction of normal incidence
%  beta is the polar angle measured from the y-axis
%  alpha is the azimuthal angle measured around the z-axis
spherAlt = @(alpha,beta) [cos(beta)*sin(alpha);-cos(beta)*cos(alpha);sin(beta)]; 

%% initialize variables
alphas = deg2rad(alphas);
betas = deg2rad(betas);
TavgBare = [];
TavgRef = [];
TavgCloak = [];

%% calculate angle-averaged transmittances
% no cloak, no reference structure
n = [0 1 0];
[BETAS,ALPHAS] = meshgrid(betas,alphas);
integrand = nan*ones(size(BETAS));
intNorm = nan*ones(size(BETAS));
for m = 1:size(BETAS,1)
    for k = 1:size(BETAS,2)
        if reflections
            [~,T] = Fresnel(spherAlt(ALPHAS(m,k),BETAS(m,k))',n,n1,n2,'none');
        else
            T = 1;
        end
        integrand(m,k) = cos(BETAS(m,k)) * T;
        intNorm(m,k) = cos(BETAS(m,k));
    end
end
TavgBare = trapz(alphas,trapz(betas,integrand,2)) / trapz(alphas,trapz(betas,intNorm,2)); % angle-averaged transmission with cloak (isotropic illumination)

% with cloak
n = [0 1 0];
[BETAS,ALPHAS] = meshgrid(betas,alphas);
integrand = nan*ones(size(BETAS));
intNorm = nan*ones(size(BETAS));
for m = 1:size(BETAS,1)
    for k = 1:size(BETAS,2)
        if reflections
            [~,T] = Fresnel(spherAlt(ALPHAS(m,k),BETAS(m,k))',n,n1,n2,'none');
        else
            T = 1;
        end
        integrand(m,k) = cos(BETAS(m,k)) * T * (P(rad2deg(ALPHAS(m,k)),rad2deg(BETAS(m,k))));
        intNorm(m,k) = cos(BETAS(m,k));
    end
end

TavgCloak = trapz(alphas,trapz(betas,integrand,2)) / trapz(alphas,trapz(betas,intNorm,2)); % angle-averaged transmission with cloak (isotropic illumination)


if referenced
    % with reference structure
    n = [0 1 0];
    [BETAS,ALPHAS] = meshgrid(betas,alphas);
    integrand = nan*ones(size(BETAS));
    intNorm = nan*ones(size(BETAS));
    for m = 1:size(BETAS,1)
        for k = 1:size(BETAS,2)
            if reflections
                [~,T] = Fresnel(spherAlt(ALPHAS(m,k),BETAS(m,k))',n,n1,n2,'none');
            else
                T = 1;
            end
            integrand(m,k) = cos(BETAS(m,k)) * T * Pref(rad2deg(ALPHAS(m,k)),rad2deg(BETAS(m,k)));
            intNorm(m,k) = cos(BETAS(m,k));
        end
    end

    TavgRef = trapz(alphas,trapz(betas,integrand,2)) / trapz(alphas,trapz(betas,intNorm,2)); % angle-averaged transmission with cloak (isotropic illumination)
end

