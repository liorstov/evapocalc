function [sim_vals, sim_wet] = synthetic_daily_data_wbl(L,p_wet_afterwet,p_wet_afterdry,wbl_phat)

rd = rand(L,365);
sim_wet = rand(L,365);
sim_vals = zeros(L,365);
sim_wet(:,1) = rd(:,1) < p_wet_afterdry(1);
for k = 2: 365
    sim_wet(:,k) = (sim_wet(:,k-1) & (rd(:,k) < p_wet_afterwet(k))) | ~sim_wet(:,k-1) & (rd(:,k) < p_wet_afterdry(k));
end
sim_vals(find(sim_wet)) = wblinv(rand(sum(sim_wet(:)),1), wbl_phat(1), wbl_phat(2));
    
end