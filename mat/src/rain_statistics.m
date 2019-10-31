function rain_statistics(rain_ifname,rain_fname,month_days)

% loads data
    load(rain_ifname,'R'); 
    R.name
    SMEV_1type_phat = R.SMEV_1ty.phat;

% available years
    yrs = 1948:2018;
    for y = 1: numel(yrs); tmp_annual_precip(y) = nansum(R.vals(R.time>=datenum([yrs(y) 9 1 0 0 0]) & R.time<datenum([yrs(y)+1 9 1 0 0 0]))); end
    yrs(tmp_annual_precip==0)=[]; clear tmp_annual_precip
    M = numel(yrs);
    
% type (implies wet/dry)
    obs_type = zeros(numel(yrs),365);
    obs_rain = zeros(numel(yrs),365);
    for y = 1: numel(yrs)
        k = 0; % days counter
        for mth = [10: 12, 1: 9]
            if mth>=9; yr = yrs(y); else yr = yrs(y)+1; end
            for dy = 1: month_days(mth)
                k = k + 1;
                if any(R.time==datenum([yr mth dy 0 0 0]));
                    obs_type(y,k) = R.type(R.time==datenum([yr mth dy 0 0 0]));
                    obs_rain(y,k) = R.vals(R.time==datenum([yr mth dy 0 0 0]));
                end
            end
        end
    end
    
% wet/dry probability
    obs_wet = obs_type~=0;
    %all days following a wet day
    pet_is_wet_prev = [false(M,1) obs_wet(:,1:end-1)];
    obs_p_wet = sum(obs_wet,1)/M;
    obs_p_wet_afterdry = sum(obs_wet & ~pet_is_wet_prev,1)./sum(~pet_is_wet_prev,1); obs_p_wet_afterdry(isnan(obs_p_wet_afterdry)) = 0;
    obs_p_wet_afterwet = sum(obs_wet & pet_is_wet_prev,1)./sum(pet_is_wet_prev,1); obs_p_wet_afterwet(isnan(obs_p_wet_afterwet)) = 0;

save(rain_ifname,'R','obs_rain','obs_type','obs_wet','obs_p_wet','obs_p_wet_afterdry','obs_p_wet_afterwet','yrs','M','SMEV_1type_phat');

end