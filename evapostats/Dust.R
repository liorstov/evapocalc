



DustTrop = read_csv("dust.csv", na = "n.d.") %>% dplyr::select(1:6, Ca) %>%  mutate(deployment = dmy_hm(deployment), recovery = dmy_hm(recovery), ppm = Ca/DustLoad*1000/10000) 


Dustsum = DustTrop  %>% group_by(Fraction)%>%  summarise( mean(ppm), std(ppm)) 

sum(DustTrop[1,7:29])/1000

