classdef irradianceDataClass < handle
    properties
        inputFolder
        monthlyData
        dailyData
        interpolatedData
        interpolatedDataDiff
        fromDirectNormalIrradiance;
    end
    
    methods
        function obj = irradianceDataClass(input_folder,fromDNI)
            obj.inputFolder = input_folder;
            obj.fromDirectNormalIrradiance = fromDNI;

            fid = fopen([input_folder filesep 'monthly.txt']);
            tline = fgetl(fid);
            takeData = false; m = 1;
            while ischar(tline)
                if takeData && m<=12
                    C = strsplit(tline,'\t');
                    obj.monthlyData(m).Hh = str2double(C{2});
                    obj.monthlyData(m).DG = str2double(C{3});
                    m = m+1;
                elseif ~isempty(regexp(tline,'Month\t+Hh\t+D/G', 'once'))
                    takeData = true;
                end             

                tline = fgetl(fid);
            end
            
            fclose(fid);
                       
            obj.dailyData = {};
            ird = []; irdDiff = [];
            P = [];
            for month = 1:12
                dd = dailyIrradianceDataClass(input_folder,month);
                obj.dailyData{month} = dd;
                
                P = [P; repmat(date2doy(15,month),length(dd.localSunTimes),1), dd.localSunTimes'];
                if obj.fromDirectNormalIrradiance
                    ird = [ird; dd.directNormalIrradiance(dd.localSunTimes)'];  % from direct normal irradiance
                else
                    ird = [ird; dd.globalIrradiance(dd.localSunTimes)'-dd.diffuseIrradiance(dd.localSunTimes)']; % from global irradiance
                    irdDiff = [irdDiff; dd.diffuseIrradiance(dd.localSunTimes)'];
                end                
            end
            obj.interpolatedData = scatteredInterpolant(P, ird, 'linear','nearest');
            if ~obj.fromDirectNormalIrradiance
                obj.interpolatedDataDiff = scatteredInterpolant(P, irdDiff, 'linear', 'nearest');
            end
        end
    end
    
end