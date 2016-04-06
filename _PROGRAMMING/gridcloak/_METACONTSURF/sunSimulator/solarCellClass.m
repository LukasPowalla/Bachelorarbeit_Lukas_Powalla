classdef solarCellClass < handle
    properties (Constant, Access = private)
        ignoreSunHeightsBelow = 0
    end
    properties
        azimuth = 0     % cell azimuth
                        %  0 = cell surface faces south (for positive inclination)
                        %  90 = cell surface faces east (for positive
                        %       inclination)
        inclination = 0 % cell inclination
        rotation = 0    % cell rotation on roof around surface normal of the cell.
                        %  0 = wires going along the slope of the roof
    end
   
    methods (Static, Access = private)
        function result = R(phi,u)
        % rotation matrix for active rotation by angle phi around axis u
           u = u/norm(u);
           result = [...
            cosd(phi)+u(1)^2*(1-cosd(phi)), u(1)*u(2)*(1-cosd(phi))-u(3)*sind(phi), u(1)*u(3)*(1-cosd(phi))+u(2)*sind(phi);
            u(2)*u(1)*(1-cosd(phi))+u(3)*sind(phi), cosd(phi)+u(2)^2*(1-cosd(phi)), u(2)*u(3)*(1-cosd(phi))-u(1)*sind(phi);
            u(3)*u(1)*(1-cosd(phi))-u(2)*sind(phi), u(3)*u(2)*(1-cosd(phi))+u(1)*sind(phi), cosd(phi)+u(3)^2*(1-cosd(phi))]; 
        end
    end
    
    methods
        function obj = solarCellClass(azimuth, inclination, rotation)
            obj.azimuth = azimuth;
            obj.inclination = inclination;
            obj.rotation = rotation;
        end
        
        function [alpha,beta] = incidentAngles(obj, city, dayOfYear, localTime)
            % returns incident angles relative to solar cell orientation
            %  alpha: angle between incident light and solar cell surface normal 
            %        measured in the plane created by a vector within the 
            %        cell surface and perpendicular to the wire axis, and
            %        the cell surface normal
            %  beta: angle between incident light and solar cell surface normal 
            %        measured in the plane created by the wire axis and the 
            %        cell surface normal
            [sunHeight,sunAzimuth] = city.sunAngles(dayOfYear,localTime);
            cellInclination = obj.inclination;
            cellAzimuth = obj.azimuth; 
            cellRotation = obj.rotation;
            R = @solarCellClass.R;
            
            v2 = R(cellAzimuth,[0;0;1])*[0;1;0];
            n = R(-cellInclination,v2)*[0;0;1];
            v1 = R(cellRotation,n)*R(-cellInclination,v2)*R(cellAzimuth,[0;0;1])*[1;0;0];
            v2 = R(cellRotation,n)*v2;

            s = [cosd(sunAzimuth).*cosd(sunHeight); -sind(sunAzimuth).*cosd(sunHeight); sind(sunHeight)];
 
            
            v1s = repmat(v1,1,size(s,2));
            v2s = repmat(v2,1,size(s,2));
            ns = repmat(n,1,size(s,2));
            
            v1projs = s-repmat(dot(s,v1s),3,1).*v1s;
            v1projs = v1projs ./ repmat(sqrt(sum(v1projs.^2)),3,1);
            
            v2projs = s-repmat(dot(s,v2s),3,1).*v2s;
            v2projs = v2projs ./ repmat(sqrt(sum(v2projs.^2)),3,1);
            
            alpha = acosd(dot(v1projs,ns));
            beta = acosd(dot(v2projs,ns));
        end
        
        function [result,theta] = inclinationFactor(obj, city, dayOfYear, localTime)
            [sunHeight,sunAzimuth] = city.sunAngles(dayOfYear,localTime);
            cellInclination = obj.inclination;
            cellAzimuth = obj.azimuth; 
            theta = acosd(sind(cellInclination)*(-cosd(sunHeight)).*cosd(sunAzimuth-cellAzimuth)+sind(sunHeight)*cosd(cellInclination));
            if city.irradianceData.fromDirectNormalIrradiance
                result = cosd(theta);
            else
                result = cosd(theta)./sind(sunHeight);   
            end
            result(theta > 90 | theta < -90 | sunHeight < solarCellClass.ignoreSunHeightsBelow) = 0;
        end
        
        function [result,resultDiff,theta,inclinationFactor] = totalIntensity(obj, city, dayOfYear, localTime)
             result = []; resultDiff = [];
             [inclinationFactor,theta] = obj.inclinationFactor(city,dayOfYear,localTime);
             ir = city.irradianceData.interpolatedData;
             irDiff = city.irradianceData.interpolatedDataDiff;
             if ~isempty(ir)
                result = ir(repmat(dayOfYear,1,length(localTime)),city.localSunTime(dayOfYear,localTime)).*inclinationFactor;
             end
             if ~isempty(irDiff)
                resultDiff = irDiff(repmat(dayOfYear,1,length(localTime)),city.localSunTime(dayOfYear,localTime)).*0.5*(1+cosd(obj.inclination));
             end
        end
    end
    
    
end