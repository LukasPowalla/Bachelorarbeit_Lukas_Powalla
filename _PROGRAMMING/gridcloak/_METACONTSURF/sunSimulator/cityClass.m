classdef cityClass < handle
    properties
        description = ''
        timeZone = 0
        latitude = 0
        longitude = 0
        irradianceData
    end
    methods (Static)
        function result = declination(dayOfYear)
            J = 360*dayOfYear/365;
            result = 0.3948-23.2559*cosd(J+9.1)-0.3915*cosd(2*J+5.4)-0.1764*cos(2*J+26);  % declination in degrees
        end
    end
    methods
        function obj = cityClass(varargin)
            % constructor
            %  cityClass(description)
            %  cityClass(timeZone, latitude, longitude)
            %  cityClass(description, timeZone, latitude, longitude)
            %  cityClass(description, timeZone, latitude, longitude, irradianceDataInputFolder, irradianceDataFromDNI)
            switch nargin
                case 1
                    obj.description = varargin{1};
                case 3
                    obj.timeZone = varargin{1};
                    obj.latitude = varargin{2};
                    obj.longitude = varargin{3};
                case 4
                    obj.description = varargin{1};
                    obj.timeZone = varargin{2};
                    obj.latitude = varargin{3};
                    obj.longitude = varargin{4};
                case 6
                    obj.description = varargin{1};
                    obj.timeZone = varargin{2};
                    obj.latitude = varargin{3};
                    obj.longitude = varargin{4};
                    obj.irradianceData = irradianceDataClass(varargin{5},varargin{6});
                otherwise
                    error('Either 1, 3, 4, or 6 arguments expected.');
            end
        end

        function result = localSunTime(obj,dayOfYear,localTime)
            J = 360*dayOfYear/365;
            timeEq = (0.0066+7.3525*cosd(J+85.9)+9.9359*cosd(2*J+108.9)+0.3387*cosd(3*J+105.2))/60; % time equation in hours
            meanLocalTime = localTime-obj.timeZone+4*obj.longitude/180*pi;     % MOZ in hours
            result = meanLocalTime + timeEq;             % WOZ in hours
        end

        function [sunHeight, sunAzimuth] = sunAngles(obj,dayOfYear,localTime)
            trueLocalTime = obj.localSunTime(dayOfYear,localTime);
            declination = cityClass.declination(dayOfYear);
            
            hourAngle = (12-trueLocalTime)*15;
            sunHeight = asind(cosd(hourAngle)*cosd(obj.latitude)*cosd(declination)+sind(obj.latitude)*sind(declination));
            sunAzimuth_ = acosd((sind(sunHeight)*sind(obj.latitude)-sind(declination))./(cosd(sunHeight)*cosd(obj.latitude)));
            
            sunAzimuth = 180 + sunAzimuth_;
            sunAzimuth(trueLocalTime <= 12 & trueLocalTime > 0) = 180 - sunAzimuth_(trueLocalTime <= 12 & trueLocalTime > 0);
            sunAzimuth(trueLocalTime <= 0) = sunAzimuth_(trueLocalTime <= 0) - 180;
        end
        
    end
end