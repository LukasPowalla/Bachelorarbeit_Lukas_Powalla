classdef dailyIrradianceDataClass < handle
    properties (Access = private)
       data 
    end
    
    properties (SetAccess = private)
        localSunTimes
        inputFolder
        month
        globalIrradiance
        diffuseIrradiance
        clearskyIrradiance
        directNormalIrradiance
        clearskyNormalIrradiance    
    end
    
    methods
        function obj = dailyIrradianceDataClass(input_folder,month)
            obj.inputFolder = input_folder;
            obj.month = month;
            input_file = [input_folder filesep sprintf('daily_%02.0f.txt',month)];
            
            d = importdata(input_file,'\t');
            
            obj.data = []; n = 1;
            for m = 1:size(d.data,2)
                col = d.data(:,m);
                if ~max(isnan(col))
                    obj.data(:,n) = col;
                    n = n+1;
                end
            end
           
            fid = fopen(input_file,'r');
            tline = fgetl(fid);
            m = 1; obj.localSunTimes = [];
            while ischar(tline)
                r = regexp(tline,'[0-9]{2}\:[0-9]{2}','match');
                if ~isempty(r)
                    C = strsplit(r{1},':');
                    hh = str2double(C{1}); mm = str2double(C{2});
                    obj.localSunTimes(m) = hh+mm/60;
                    m = m+1;
                end
                tline = fgetl(fid);
            end
            fclose(fid);
            
            dt = min(diff(obj.localSunTimes));
            t1 = obj.localSunTimes(1); tend = obj.localSunTimes(end);
            obj.localSunTimes = [t1-dt obj.localSunTimes tend+dt];
            
            obj.data = [zeros(1,size(obj.data,2)); obj.data; zeros(1,size(obj.data,2))];
            
            obj.globalIrradiance = griddedInterpolant(obj.localSunTimes,obj.data(:,1),'linear','nearest');
            obj.diffuseIrradiance = griddedInterpolant(obj.localSunTimes,obj.data(:,2),'linear','nearest');
            obj.clearskyIrradiance = griddedInterpolant(obj.localSunTimes,obj.data(:,3),'linear','nearest');
            obj.directNormalIrradiance = griddedInterpolant(obj.localSunTimes,obj.data(:,4),'linear','nearest');
            obj.clearskyNormalIrradiance = griddedInterpolant(obj.localSunTimes,obj.data(:,5),'linear','nearest'); 
        end
        
    end

    
end
