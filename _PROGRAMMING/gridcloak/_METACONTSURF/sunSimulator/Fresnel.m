function [R,T] = Fresnel(v,n,n1,n2,pol)
% calculates the Fresnel Reflection and Transmission of a ray with velocity
% v approaching a surface with surface normal n, coming from a medium with
% refractive index n1, and transmitting into a medium with refractive index
% n2. Its polarization is taken to be pol (can be 's' or 'p', or anything
% else for unpolarized light).
% v has to be row vector or a matrix whose rows are the individual v
% vectors. n has to be a row vector, n1 and n2 must be scalars.

    switch size(v,2)
        case 2
            norm_v = (v(:,1).^2+v(:,2).^2).^0.5;
        case 3
            norm_v = (v(:,1).^2+v(:,2).^2+v(:,3).^2).^0.5;
        otherwise
            error('v must be a row vector or matrix of row vectors with either 2 or 3 columns');
    end
    phi1 = acos(abs(v*n')./norm_v/norm(n));
    
    phi2 = asin(n1/n2*sin(phi1));
    
    rs = (n1*cos(phi1)-n2*cos(phi2))./(n1*cos(phi1)+n2*cos(phi2));
    rp = (n2*cos(phi1)-n1*cos(phi2))./(n1*cos(phi2)+n2*cos(phi1));
    
    switch pol
        case 's'
            R = abs(rs).^2;
        case 'p'
            R = abs(rp).^2;
        otherwise
            R = (abs(rs).^2 + abs(rp).^2) * .5;
    end
    
    T = 1-R;


