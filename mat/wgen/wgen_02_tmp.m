% Point rainfall generator conditioned on seasonality and synoptic type
% some parts based on the SMEV analyses
% 
% Notes:
%   *   feb 29 in leap years are ignored
% 
% % 
% 
%   F Marra, Sep 2019
%
% 
%%

clear
dbpath = 'C:\Users\marraf\Dropbox\';
wpath = [dbpath,'work\israel\HUJI_CC_proj\weather_gen\'];
synoptics = [dbpath,'work\papers_nostri\MEV_CC\analyses\synoptic\synoptic_classification_1948-10_May_2018_NewClasses_Lows_RST-Other_AfterAggregation']; % new classification for paper 2
month_days = [31 28 31 30 31 30 31 31 30 31 30 31];
month_str = {'Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug'};
TYPES_col   = {[0 .5 1],[.8 .2 0],'k'};
sz = 30;
sz2 = 15;
% id = 244730; nm = 'jer'; % best choice for RAIN and PET data availability
id = 347700; nm = 'elat';
L_multip = 50; % n. of simulated series of M years
% L = 69*1000; % n. of simulated years

%% loads daily archive for general climatic info
fname = [wpath,nm,'_data'];
if exist([fname,'.mat'],'file')~=2
    load([dbpath,'work\papers_nostri\MEV_CC\analyses\05_35yrs_paper2\relevant_35_elabs_1000.mat']);
    R = R([R.id]==id); R.type(R.type==0) = 2; save(fname,'R'); R.name
else
    load(fname,'R'); 
    R.name
end

%% analysis of observations

% available years
    yrs = 1948:2017;
    for y = 1: numel(yrs); tmp_annual_precip(y) = nansum(R.vals(R.time>=datenum([yrs(y) 9 1 0 0 0]) & R.time<datenum([yrs(y)+1 9 1 0 0 0]))); end
    yrs(tmp_annual_precip==0)=[]; clear tmp_annual_precip
    M = numel(yrs);
    
% type (implies wet/dry)
    obs_type = zeros(numel(yrs),365);
    obs_vals = zeros(numel(yrs),365);
    for y = 1: numel(yrs)
        k = 0; % days counter
        for mth = [9: 12, 1: 8]
            if mth>=9; yr = yrs(y); else yr = yrs(y)+1; end
            for dy = 1: month_days(mth)
                k = k + 1;
                if any(R.time==datenum([yr mth dy 0 0 0]));
                    obs_type(y,k) = R.type(R.time==datenum([yr mth dy 0 0 0]));
                    obs_vals(y,k) = R.vals(R.time==datenum([yr mth dy 0 0 0]));
                end
            end
        end
    end
    
% wet/dry probability
    obs_wet = obs_type~=0;
    obs_wet_prev = [false(M,1) obs_wet(:,1:end-1)];
    obs_p_wet = sum(obs_wet,1)/M;
    obs_p_wet_afterdry = sum(obs_wet & ~obs_wet_prev,1)./sum(~obs_wet_prev,1); obs_p_wet_afterdry(isnan(obs_p_wet_afterdry)) = 0;
    obs_p_wet_afterwet = sum(obs_wet & obs_wet_prev,1)./sum(obs_wet_prev,1); obs_p_wet_afterwet(isnan(obs_p_wet_afterwet)) = 0;    
    
%% weather generator

% wet/dry probability for simulation (smoothing)
    p_wet = conv(conv(obs_p_wet,ones(1,sz)/sz,'same'),ones(1,sz2)/sz2,'same');
    p_wet_afterdry = conv(conv(obs_p_wet_afterdry,ones(1,sz)/sz,'same'),ones(1,sz2)/sz2,'same');
    p_wet_afterwet = conv(conv(obs_p_wet_afterwet,ones(1,sz)/sz,'same'),ones(1,sz2)/sz2,'same');

% generates synth data
    tic
    L = M*L_multip;
%     L = 100000;
    [sim_vals, sim_wet] = synthetic_daily_data_wbl(L,p_wet_afterwet,p_wet_afterdry,R.SMEV_1ty.phat);
    toc

%% check simulation statistics

[obs_annual_precip,obs_ams,obs_annual_n,obs_n_events,obs_durations] = calc_annual_stats(obs_wet,obs_vals);
[sim_annual_precip,sim_ams,sim_annual_n,sim_n_events,sim_durations] = calc_annual_stats(sim_wet,sim_vals);
sim_p_wet = sum(sim_wet,1)/L;

% distrib stats for uncertainty
tic
sim_annual_precip_unc = nan(L_multip,M);
sim_ams_unc = nan(L_multip,M);
sim_annual_n_unc = nan(L_multip,M);
sim_n_events_unc = nan(L_multip,M);
sim_durations_unc = cell(L_multip,1);
sim_p_wet_unc = nan(L_multip,365);
for l = 1: L_multip
    [sim_annual_precip_unc(l,:),sim_ams_unc(l,:),sim_annual_n_unc(l,:),sim_n_events_unc(l,:),sim_durations_unc{l}] = ...
        calc_annual_stats(sim_wet( (l-1)*M + (1:M),:), sim_vals( (l-1)*M + (1:M),:) );
    sim_p_wet_unc(l,:) = sum(sim_wet( (l-1)*M + (1:M),:))/M;
%     sim_durations_unc(l,1:numel(sdu)) = sdu;
end
toc

%% figure to check simulation statistics
dxr = 30;

ff=figure; set(ff,'position',[100 100 1200 600]);

    subplot(2,3,1);
        ll(1) = plot(obs_p_wet,'r-'); hold on
        ll(2) = plot(sim_p_wet,'k-','linewidth',2); hold on
        jbfill(1:365,quantile(sim_p_wet_unc,.95,1),quantile(sim_p_wet_unc,.05,1),'k','none',1,.2);
        ll(3) = plot(p_wet,'g--','linewidth',2); hold on
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
        
% print([fname,'_simul.png'],'-dpng','-r300'); close(ff)
