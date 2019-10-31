function synth_data = synthetic_timeseries(sim_rain,sim_pet,start_year,L)

synth_data.time = (datenum([start_year 9 1 0 0 0]): 1: datenum([start_year+L 8 31 0 0 0]));

synth_data.P = single(nan(size(synth_data.time)));
synth_data.PET = single(nan(size(synth_data.time)));
last=0;
for y = 1: L
    if leapyear(start_year+y)
        synth_data.P(last+1:last+59) = sim_rain(y,1:59);
        synth_data.P(last+60) = sim_rain(y,59);
        synth_data.P(last+61:last+366) = sim_rain(y,60:365);
        
        synth_data.PET(last+1:last+59) = sim_pet(y,1:59);
        synth_data.PET(last+60) = sim_pet(y,59);
        synth_data.PET(last+61:last+366) = sim_pet(y,60:365);
        last = last+366;
    else
        synth_data.P(last+1:last+365) = sim_rain(y,:);
        synth_data.PET(last+1:last+365) = sim_pet(y,:);
        last = last+365;
    end
    
end

end
