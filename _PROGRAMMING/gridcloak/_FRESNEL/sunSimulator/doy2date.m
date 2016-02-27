function [d,m] = doy2date(doy)
    [~,m,d,~,~,~] = datevec(doy+366);
end