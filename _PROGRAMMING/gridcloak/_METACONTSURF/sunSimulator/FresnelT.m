function T = FresnelT(angle,n,pol)
% calculates the combined Fresnel transmission through a layer stack
% consisting of layers with refractive indices specified in vector n.
% angle specifies the incidence angle and might be a vector of angles.
% pol specifies the polarization of the rays (may be 's', 'p', or anything
% else for unpolarized light. The calculated transmission T is only
% calculated from the reflectivity amplitudes, thus interference is NOT considered.
% Also, only a single transmitted ray is considered (i.e. no interior
% reflections).

    if length(n)<2
        error('n must be a vector of length 2 or more');
    end
    num_interfaces = length(n)-1;
    
    v = 1/n(1)*[sind(angle) -cosd(angle)];
    
    T = zeros(num_interfaces,length(angle));
    
    beta = angle;
    for i = 1:num_interfaces
       
        n1 = n(i);
        n2 = n(i+1);
        
        [~,T(i,:)] = Fresnel(v,[0 1],n1,n2,pol);
        beta = asind(n1/n2 * sind(beta));
        v = 1/n2*[sind(beta) -cosd(beta)];

    end
    
    T = prod(T,1);

end