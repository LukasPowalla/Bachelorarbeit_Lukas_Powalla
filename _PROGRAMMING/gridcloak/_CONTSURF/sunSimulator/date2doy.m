function result = date2doy(day,month)
    result = datenum(ones(1,length(day)),month,day)-366;
end