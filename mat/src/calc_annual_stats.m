function [annual_precip,ams,annual_n,events_n,events_durations] = calc_annual_stats(iswet,vals)
    
    annual_precip = sum(vals,2);
    ams = max(vals,[],2);
    annual_n = sum(iswet,2);

    if nargout>3
        
        wet_array = reshape(iswet',numel(iswet),1);
        yfrom = find(wet_array(2:end) & ~wet_array(1:end-1)) + 1; if wet_array(1); yfrom = cat(1,1,yfrom); end
        yto = find(~wet_array(2:end) & wet_array(1:end-1)); if wet_array(end); yto = cat(1,yto,numel(wet_array)); end
        events_durations = yto - yfrom + 1;
        
        if nargout>4
            
            isfrom = false(size(iswet)); isfrom(yfrom)=true;
            events_n = sum(isfrom,2);
        
        end
        
    end
    
end