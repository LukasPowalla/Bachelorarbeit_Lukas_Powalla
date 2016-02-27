function vtrans = Snell(v,n,n1,n2)
    %% refract ray with velocity v at plane with normal vector n separating regions with refractive indices n1 and n2
    
    % normalize input vectors and ensure n is pointing upwards (against
    % direction of v)
    v_ = repmat(sqrt(dot(v,v,2)),1,3);
    v = v./v_;
    n = -n./repmat(sqrt(dot(n,n,2)),1,3).*repmat(sign(dot(n,v,2)),1,3);
    
    % calculate direction of refracted ray
    cnv = cross(n,v,2);
    
    vtrans = nan * ones(size(v));
    
    transmit = n2/n1>sqrt(dot(cnv,cnv,2));
%     if n2/n1>sqrt(dot(cnv,cnv,2)) 
        vtrans(transmit,:) = n1/n2*(cross(n(transmit,:),cross(-n(transmit,:),v(transmit,:),2),2))-n(transmit,:).*sqrt(1-(n1/n2)^2*(repmat(dot(cnv(transmit,:),cnv(transmit,:),2),1,3)));
        % adjust velocity
        vtrans(transmit,:) = v_(transmit,:)*n1/n2.*vtrans(transmit,:);
%     else % total reflection
        if ~all(transmit)
            vtrans(~transmit,:) = 2*repmat(dot(v(~transmit,:),n(~transmit,:),2),1,3).*n(~transmit,:)+v(~transmit,:);
        end
%     end    

end