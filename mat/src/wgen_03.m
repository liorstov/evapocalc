% Point rainfall generator conditioned on seasonality and synoptic type
% some parts based on the SMEV analyses
% Note: in the definition of the statistics, Feb 29 are ignored
%       Feb 29 are produced, when needed, as replicas of Feb 28
% 
%   F Marra, Sep 2019
%
% 
clear
month_days = [31 28 31 30 31 30 31 31 30 31 30 31];
month_str = {'Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug'};

% study case
    id = 244730; nm = 'elat'; % best choice for RAIN and PET data availability
    L_multip = 1000; % L = M*L_multip; % multiplicator of the n. of simulated years
    start_year = 1999;
    
% RAIN
    rain_full_archive = ['work\papers_nostri\MEV_CC\analyses\05_35yrs_paper2\relevant_35_elabs_1000.mat'];
    rain_ifname = [nm,'_data.mat'];
    rain_stat_fname = [nm,'_rain.mat'];
    sz = 40;  %[days] % smoothing window size 1
    sz2 = 30; %[days] % smoothing window size 2

% PET
    pet_ifname = 'Elat_PET.txt';
    pet_stat_fname = [nm,'_PET.mat'];
    pet_sz = 30; %[days] 

%% preprocessing: if needed creates rain and PET statistics
% in the definition of the statistics, Feb 29 are ignored
    if exist(rain_stat_fname,'file')~=2
        if exist(rain_ifname,'file')~=2; load(rain_full_archive); R = R([R.id]==id); R.type(R.type==0) = 2; save(rain_ifname,'R'); end
        rain_statistics(rain_ifname,rain_stat_fname,month_days);
    end

    if exist(pet_stat_fname,'file')~=2
        load(rain_stat_fname,'obs_wet','yrs');
        PET_statistics(pet_ifname,pet_stat_fname,pet_sz,month_days,obs_wet,yrs);
    end
    
%% weather generator  

% loads statistics of rain and PET
    load(rain_ifname);
    load(pet_stat_fname);
    
% wet/dry probability for simulation (smoothing)
    p_wet = conv(conv(obs_p_wet,ones(1,sz)/sz,'same'),ones(1,sz2)/sz2,'same');
    p_wet_afterdry = conv(conv(obs_p_wet_afterdry,ones(1,sz)/sz,'same'),ones(1,sz2)/sz2,'same');
    p_wet_afterwet = conv(conv(obs_p_wet_afterwet,ones(1,sz)/sz,'same'),ones(1,sz2)/sz2,'same');

% generates synth rainfall data
    fprintf('generating RAIN: '); tic
        if exist('L_multip','var')==1; L = M*L_multip; fprintf('(using L_multip=%i)',L_multip); end
        [sim_rain, sim_wet] = synthetic_daily_data_wbl(L,p_wet_afterwet,p_wet_afterdry,SMEV_1type_phat);
    toc

% generates synth PET data based on input rain data
    fprintf('generating PET: '); tic
        sim_pet = synthetic_PET_data(sim_wet,pet_wet_avg,pet_wet_std,pet_dad_avg,pet_dad_std,pet_daw_avg,pet_daw_std);
    toc
    
% % saves synth data
%     save([wpath,'synthetic_rain_PET_',num2str(L),'years'],'sim_rain','sim_pet');
    
% creates realistic timeseries (also including Feb 29 in leap years as
% replicas of Feb 28)
    fprintf('building timeseries: '); tic
        synth_data = synthetic_timeseries(sim_rain,sim_pet,start_year,L);
    toc
    save(['synth_data__cwl__',num2str(L),'years'],'synth_data');
    
%% checks PET simulation statistics
% 
% load([wpath,'synthetic_rain_PET_',num2str(L),'years']);
% sim_pet(sim_pet<.5) = sim_pet(sim_pet<.5)*1.5;
% sim_pet(sim_pet<0) = -sim_pet(sim_pet<0);
% sim_pet(sim_pet<0) = NaN;

    fprintf('checking PET: '); tic
        xpet=0:.1:14;
        yobs = hist(obs_pet(:),xpet)/pet_M;
        ysim = hist(sim_pet(:),xpet)/L;
        ysim_unc = nan(fix(L/pet_M),numel(xpet));
        for l = 1: fix(L/pet_M)
            tmp = sim_pet(((l-1)*pet_M)+1:l*pet_M,:);
            ysim_unc(l,:) = hist(tmp(:),xpet)/pet_M;
        end
    toc
    
    ff=figure; set(ff,'position',[100 100 600 900]);
        subplot(2,1,1);
            plot(xpet,yobs,'r-'); hold on
            plot(xpet,ysim,'k-','linewidth',2);
            jbfill(xpet,quantile(ysim_unc,.95,1),quantile(ysim_unc,.05,1),'k','none',1,.2);
            legend({'Observed PET','Simulated PET'});
            xlim([0 14]); xlabel('PET [mm day^{-1}]'); ylabel('# of yearly occurrences'); box on
        subplot(2,1,2);
            q1 = .05;
            plot(1:365,nanmean(obs_pet,1),'r-','linewidth',2); hold on
            jbfill(1:365,quantile(obs_pet,1-q1,1),quantile(obs_pet,q1,1),'r','none',1,.2);
            plot(1:365,nanmean(sim_pet,1),'k-','linewidth',2);
            jbfill(1:365,quantile(sim_pet,1-q1,1),quantile(sim_pet,q1,1),'k','none',1,.2);
%             jbfill(1:365,quantile(sim_pet,.99,1),quantile(sim_pet,.01,1),'k','none',1,.2);
            mloc = [1 1+cumsum(month_days(1:11))];
            for mth = 1:12; text(mloc(mth)+5,-.01,month_str{mth},'rotation',45,'fontsize',10,'horizontalalignment','right'); end
            xlim([1 365]); set(gca,'xtick',mloc,'xticklabel',''); ylabel('PET [mm day^{-1}]'); box on
    print([nm,'_PET_simul2.png'],'-dpng','-r300'); 
    
%% check RAIN simulation statistics

[obs_annual_precip,obs_ams,obs_annual_n,obs_n_events,obs_durations] = calc_annual_stats(obs_wet,obs_rain);
[sim_annual_precip,sim_ams,sim_annual_n,sim_n_events,sim_durations] = calc_annual_stats(sim_wet,sim_rain);
sim_p_wet = sum(sim_wet,1)/L;

% distrib stats for uncertainty
fprintf('checking RAIN: ');tic
    sim_annual_precip_unc = nan(L_multip,M);
    sim_ams_unc = nan(L_multip,M);
    sim_annual_n_unc = nan(L_multip,M);
    sim_n_events_unc = nan(L_multip,M);
    sim_durations_unc = cell(L_multip,1);
    sim_p_wet_unc = nan(L_multip,365);
    for l = 1: L_multip
        [sim_annual_precip_unc(l,:),sim_ams_unc(l,:),sim_annual_n_unc(l,:),sim_n_events_unc(l,:),sim_durations_unc{l}] = ...
            calc_annual_stats(sim_wet( (l-1)*M + (1:M),:), sim_rain( (l-1)*M + (1:M),:) );
        sim_p_wet_unc(l,:) = sum(sim_wet( (l-1)*M + (1:M),:))/M;
    %     sim_durations_unc(l,1:numel(sdu)) = sdu;
    end
toc

dxr = 30;

ff=figure; set(ff,'position',[100 100 1200 600]);

    subplot(2,3,1);
        ll(1) = plot(obs_p_wet,'r-','linewidth',2); hold on
        ll(2) = plot(sim_p_wet,'k-','linewidth',2); hold on
        jbfill(1:365,quantile(sim_p_wet_unc,.95,1),quantile(sim_p_wet_unc,.05,1),'k','none',1,.2);
        %ll(3) = plot(p_wet,'g--','linewidth',2); hold on
        mloc = [1 1+cumsum(month_days(1:11))];
        for mth = 1:12; text(mloc(mth)+5,-.01,month_str{mth},'rotation',45,'fontsize',10,'horizontalalignment','right'); end
        legend(ll,{'Observed','Simulated','Planned'}); clear ll
        title([R.name,'; M=',num2str(M),'; L=',num2str(L)]); ylabel('Wet probability p [-]'); xlim([1 365]);
        set(gca,'xtick',mloc,'xticklabel','');
        box on
        
    subplot(2,3,2);
        map_hist_x=0:max(obs_annual_precip*1.5)/dxr:max(obs_annual_precip*1.5);
        obs_ams_hist=hist(obs_annual_precip,map_hist_x)/M;
        sim_ams_hist=hist(sim_annual_precip,map_hist_x)/L;
        sim_ams_hist_unc = nan(L_multip,numel(map_hist_x));
        for l = 1: L_multip; sim_ams_hist_unc(l,:)=hist(sim_annual_precip_unc(l,:),map_hist_x)/numel(sim_annual_precip_unc(l,:)); end
        jbfill(map_hist_x,quantile(sim_ams_hist_unc,.95,1),quantile(sim_ams_hist_unc,.05,1),'k','none',1,.2); hold on
        ll(1) = plot(map_hist_x,obs_ams_hist,'r-','linewidth',2); hold on
        ll(2) = plot(map_hist_x,sim_ams_hist,'k-','linewidth',2);
        title('Annual precipitation'); ylabel('pdf'); xlabel('Annual precipitation [mm year^{-1}]'); set(gca,'ytick',[]); xlim([map_hist_x(1) map_hist_x(end)]);
        legend(ll,{'Observed','Simulated'}); clear ll
        box on
        
    subplot(2,3,3);
        max_hist_x=0:10:max(obs_ams*1.5);
%         max_hist_x=0:max(obs_ams*1.5)/dxr:max(obs_ams*1.5);
        obs_max_hist=hist(obs_ams,max_hist_x)/M;
        sim_max_hist=hist(sim_ams,max_hist_x)/L;
        sim_max_hist_unc = nan(L_multip,numel(max_hist_x));
        for l = 1: L_multip; sim_max_hist_unc(l,:)=hist(sim_ams_unc(l,:),max_hist_x)/numel(sim_ams_unc(l,:)); end
        jbfill(max_hist_x,quantile(sim_max_hist_unc,.95,1),quantile(sim_max_hist_unc,.05,1),'k','none',1,.2); hold on
        ll(1) = plot(max_hist_x,obs_max_hist,'r-','linewidth',2); hold on
        ll(2) = plot(max_hist_x,sim_max_hist,'k-','linewidth',2);
        title('Annual maxima'); ylabel('pdf'); xlabel('Annual maxima [mm day^{-1}]'); set(gca,'ytick',[]); xlim([max_hist_x(1) max_hist_x(end)]);
        legend(ll,{'Observed','Simulated'}); clear ll
        box on
        
    subplot(2,3,4);
        n_hist_x=0:max(obs_annual_n*1.5)/dxr:max(obs_annual_n*1.5);
%         n_hist_x=0:10:max(obs_annual_n*1.5);
        obs_n_hist=hist(obs_annual_n,n_hist_x)/M;
        sim_n_hist=hist(sim_annual_n,n_hist_x)/L;
        sim_n_hist_unc = nan(L_multip,numel(n_hist_x));
        for l = 1: L_multip; sim_n_hist_unc(l,:)=hist(sim_annual_n_unc(l,:),n_hist_x)/numel(sim_annual_n_unc(l,:)); end
        jbfill(n_hist_x,quantile(sim_n_hist_unc,.95,1),quantile(sim_n_hist_unc,.05,1),'k','none',1,.2); hold on
        ll(1) = plot(n_hist_x,obs_n_hist,'r-','linewidth',2); hold on
        ll(2) = plot(n_hist_x,sim_n_hist,'k-','linewidth',2);
        title('# of annual wet days'); ylabel('pdf'); xlabel('# of annual wet days [-]'); set(gca,'ytick',[]); xlim([n_hist_x(1) n_hist_x(end)]);
        legend(ll,{'Observed','Simulated'}); clear ll
        box on
        
    subplot(2,3,5);
        evn_hist_x=0:1:max(obs_n_events*1.5);
        obs_evn_hist=hist(obs_n_events,evn_hist_x)/M;
        sim_evn_hist=hist(sim_n_events,evn_hist_x)/L;
        sim_evn_hist_unc = nan(L_multip,numel(evn_hist_x));
        for l = 1: L_multip; sim_evn_hist_unc(l,:)=hist(sim_n_events_unc(l,:),evn_hist_x)/numel(sim_n_events_unc(l,:)); end
        jbfill(evn_hist_x,quantile(sim_evn_hist_unc,.95,1),quantile(sim_evn_hist_unc,.05,1),'k','none',1,.2); hold on
        ll(1) = plot(evn_hist_x,obs_evn_hist,'r-','linewidth',2); hold on
        ll(2) = plot(evn_hist_x,sim_evn_hist,'k-','linewidth',2);
        title('# of annual events'); ylabel('pdf'); xlabel('# of annual events [-]'); set(gca,'ytick',[]); xlim([evn_hist_x(1) evn_hist_x(end)]);
        legend(ll,{'Observed','Simulated'}); clear ll
        box on
        
    subplot(2,3,6);
        dur_hist_x=0:1:max(obs_durations*1.2);
        obs_dur_hist=hist(obs_durations,dur_hist_x)/numel(obs_durations);
        sim_dur_hist=hist(sim_durations,dur_hist_x)/numel(sim_durations);
        sim_dur_hist_unc = nan(L_multip,numel(dur_hist_x));
        for l = 1: L_multip; sim_dur_hist_unc(l,:)=hist(sim_durations_unc{l},dur_hist_x)/numel(sim_durations_unc{l}); end
        jbfill(dur_hist_x,quantile(sim_dur_hist_unc,.95,1),quantile(sim_dur_hist_unc,.05,1),'k','none',1,.2); hold on
        ll(1) = plot(dur_hist_x,obs_dur_hist,'r-','linewidth',2); hold on
        ll(2) = plot(dur_hist_x,sim_dur_hist,'k-','linewidth',2);
        title('Event duration'); ylabel('pdf'); xlabel('Event duration [days]'); set(gca,'ytick',[]); xlim([dur_hist_x(1) dur_hist_x(end)]);
        legend(ll,{'Observed','Simulated'}); clear ll
        box on
        
print([nm,'_rain_simul.png'],'-dpng','-r300'); 
