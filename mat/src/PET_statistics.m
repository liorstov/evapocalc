function ofname = PET_statistics(ifname,ofname,pet_sz,month_days,obs_wet,yrs)

% loads PET data
    load('ElatPET.mat')    
    pet_time = datenum(ElatPET.time);
    pet_data = ElatPET.measure;
    pet_flag = ElatPET.isnullflag1;
    pet_data(pet_flag~=0) = NaN;
%     clear read_data fid

% available years
    pet_yrs = 1999:2018;
    for y = 1: numel(pet_yrs); tmp_n_annual_pet(y) = nansum(pet_flag~=1 & (pet_time>=datenum([pet_yrs(y) 9 1 0 0 0]) & pet_time<datenum([pet_yrs(y)+1 9 1 0 0 0]))); end
    pet_yrs(tmp_n_annual_pet<300)=[]; % clear tmp_n_annual_pet
    pet_M = numel(pet_yrs);

% observed PET
    obs_pet = nan(numel(pet_yrs),365);
    for y = 1: pet_M
        k = 0; % days counter
        for mth = [9: 12, 1: 8]
            if mth>=9; yr = pet_yrs(y); else yr = pet_yrs(y)+1; end
            for dy = 1: month_days(mth)
                k = k + 1;
                if any(pet_time==datenum([yr mth dy 0 0 0]));
                    obs_pet(y,k) = pet_data(pet_time==datenum([yr mth dy 0 0 0]));
                end
            end
        end
    end
    
% conditions: wet, dad, daw (uses info from the rain station: obs_wet)
    pet_is_wet = nan(size(obs_pet));
    for y = 1: pet_M; pet_is_wet(y,:) = obs_wet(yrs==pet_yrs(y),:); end % wet
    pet_is_wet_prev = [false(pet_M,1) pet_is_wet(:,1:end-1)];           % previous day wet
    pet_is_dad = ~pet_is_wet & ~pet_is_wet_prev;                        % dad: dry after dry
    pet_is_daw = ~pet_is_wet & pet_is_wet_prev;                         % daw: dry after wet
    pet_obs_wet = obs_pet; pet_obs_wet(~pet_is_wet) = NaN;
    pet_obs_dad = obs_pet; pet_obs_dad(~pet_is_dad) = NaN;
    pet_obs_daw = obs_pet; pet_obs_daw(~pet_is_daw) = NaN;

% seasonal statistics conditioned on wet/dad/daw
    pet_wet_avg = nan(365,1);   pet_dad_avg = nan(365,1);   pet_daw_avg = nan(365,1);
    pet_wet_std = nan(365,1);   pet_dad_std = nan(365,1);   pet_daw_std = nan(365,1);
    pet_wet_n = nan(365,1);     pet_dad_n = nan(365,1);     pet_daw_n = nan(365,1); 
    pet_hsz = round(pet_sz/2);
    for k = 1: 365
        % time interval in which the distributions are sampled
            k_from = max(mod(k-pet_hsz,365),1);
            k_to = max(mod(k+pet_hsz,365),1);
            if k_from<k_to; ks = k_from:k_to; else ks = [k_from:365 1:k_to]; end
        % wet
            pet_wet_vals = pet_obs_wet(:,ks); pet_wet_vals = pet_wet_vals(:); pet_wet_vals = pet_wet_vals(~isnan(pet_wet_vals));
            pet_wet_avg(k) = nanmean(pet_wet_vals); pet_wet_std(k) = nanstd(pet_wet_vals); pet_wet_n(k) = numel(pet_wet_vals);
        % dry after dry
            pet_dad_vals = pet_obs_dad(:,ks); pet_dad_vals = pet_dad_vals(:); pet_dad_vals = pet_dad_vals(~isnan(pet_dad_vals));
            pet_dad_avg(k) = nanmean(pet_dad_vals); pet_dad_std(k) = nanstd(pet_dad_vals); pet_dad_n(k) = numel(pet_dad_vals);
        % dry after wet
            pet_daw_vals = pet_obs_daw(:,ks); pet_daw_vals = pet_daw_vals(:); pet_daw_vals = pet_daw_vals(~isnan(pet_daw_vals));
            pet_daw_avg(k) = nanmean(pet_daw_vals); pet_daw_std(k) = nanstd(pet_daw_vals); pet_daw_n(k) = numel(pet_daw_vals);
        % % checks if normal distr is good approx for PET values: YES
        % xxx=0:.5:12;
        % ff=figure; set(ff,'visible','off');
        %     yyy = hist(pet_wet_vals,xxx);
        %     plot(xxx,yyy,'-'); hold on
        %     yyy = hist(pet_dad_vals,xxx);
        %     plot(xxx,yyy,'-'); hold on
        %     yyy = hist(pet_daw_vals,xxx);
        %     plot(xxx,yyy,'-'); hold on
        %     legend({'wet','dad','daw'});
        %     xlim([0 14]); ylim([0 70]);
        % print([wpath,'tmp\',num2str(k,'%03i')],'-dpng','-r100'); close(ff);
        % %
        clear pet_wet_vals pet_dad_vals pet_daw_vals
    end
save(ofname,'pet_yrs','pet_M', 'obs_pet','pet_is_wet','pet_is_dad','pet_is_daw',...
    'pet_wet_n','pet_wet_avg','pet_wet_std',...
    'pet_dad_n','pet_dad_avg','pet_dad_std',...
    'pet_daw_n','pet_daw_avg','pet_daw_std');


% figure; plot(pet_wet_n,'-'); hold on; plot(pet_dad_n,'-'); plot(pet_daw_n,'-'); legend({'wet','dad','daw'});
% figure; plot(pet_wet_avg,'-'); hold on; plot(pet_dad_avg,'-'); plot(pet_daw_avg,'-'); legend({'wet','dad','daw'});
% figure; plot(pet_wet_std,'-'); hold on; plot(pet_dad_std,'-'); plot(pet_daw_std,'-'); legend({'wet','dad','daw'});
    

end