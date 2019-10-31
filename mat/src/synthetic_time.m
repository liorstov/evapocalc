function sim_time = synthetic_time(L)

sim_time = (datenum([start_year 9 1 0 0 0]): 1: datenum([start_year+L 8 31 0 0 0]))';

end