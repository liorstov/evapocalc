function sim_pet = synthetic_PET_data(sim_wet,pet_wet_avg,pet_wet_std,pet_dad_avg,pet_dad_std,pet_daw_avg,pet_daw_std)

sim_wet = sim_wet==1;

L = size(sim_wet,1);

sim_wet_prev = [false(L,1) sim_wet(:,1:end-1)]; % previous day wet
sim_dad = ~sim_wet & ~sim_wet_prev;                 % dad: dry after dry
sim_daw = ~sim_wet & sim_wet_prev;                  % daw: dry after wet

wet_mu = repmat(pet_wet_avg',[L,1]); wet_std = repmat(pet_wet_std',[L,1]);
dad_mu = repmat(pet_dad_avg',[L,1]); dad_std = repmat(pet_dad_std',[L,1]);
daw_mu = repmat(pet_daw_avg',[L,1]); daw_std = repmat(pet_daw_std',[L,1]);

mu = nan(size(sim_wet)); mu(sim_wet) = wet_mu(sim_wet); mu(sim_dad) = dad_mu(sim_dad); mu(sim_daw) = daw_mu(sim_daw);
st = nan(size(sim_wet)); st(sim_wet) = wet_std(sim_wet); st(sim_dad) = dad_std(sim_dad); st(sim_daw) = daw_std(sim_daw);

sim_pet = single(norminv( rand(L,365), mu, st));
sim_pet(isnan(sim_pet))=0;

end 