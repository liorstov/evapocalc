require(parallel)
load("synth.RData")
cl <- makeCluster(5, type = "PSOCK")

clusterExport(cl, c("Zel11Observed", "T1.9Observed", "T1.10Observed", "T2.1Observed", "T1.1Observed", "opt.AETF", "opt.dust", "opt.FC", "opt.sulfate", "opt.WP", "Zel1Observed"))

clusterEvalQ(cl, {
    require(tidyverse)
    require(Rcpp)
    require(ggplot2)
    require(reshape2)
    require(zoo)
    require(fitdistrplus)

    require(tictoc)
    Rcpp::sourceCpp('C:/Users/liorst/source/repos/evapocalc/Calcyp/CSM.cpp', verbose = TRUE, rebuild = 0);
    b <<- new(CSMCLASS);
    source("C:/Users/liorst/source/repos/evapocalc/evapostats/Functions.R", encoding = "Windows-1252")

    opt.AETF = 1.2;
    opt.WP = 0.013;
    opt.FC = 0.1;
    opt.sulfate = 13;
    opt.dust = 6;
    seq.AETF = seq(0.8, 1.2, by = 0.4) * opt.AETF %>% rep(60);
    seq.WP = seq(0.8, 1.2, by = 0.4) * opt.WP %>% rep(60);
    seq.FC = seq(0.8, 1.2, by = 0.4) * opt.FC %>% rep(60);
    seq.rainSeq = seq(0.8, 1.2, by = 0.4) * opt.sulfate %>% rep(60);
    seq.dustSeq = seq(0.8, 1.2, by = 0.4) * opt.dust %>% rep(60);

    RainArr = seq(10, 20, by = 1);
    DustArr = c(seq(1, 26, by = 5));
    repetition = 1:5
    rainDustArray = as.matrix(crossing(RainArr, DustArr, repetition))
})
results = list()
results$T1.10 = list();
results$T1.9 = list();
results$zel11 = list();



days = seq(5,20, by = 5);
annual = seq(1,200, by = 10);
array = tibble(crossing(days, annual), mean = annual / days) %>% filter(mean >= 3 & mean <= 11)

for (i in 6:6) {
     #eilat
    IMSRain = GetImsRain(station = 347700, stationEvap = 347704);
    rainSeriesResults = GenerateSeries(NumOfSeries = 900, IMSRain = IMSRain, AnuualRain = array$annual[i], WetDays = array$days[i])

    PETresults = PETGen(rainSeriesResults$SynthRain, IMSRain, 30);
    SynthRain = rainSeriesResults$SynthRain;
    SynthRain$PET = PETresults$SynthPET;
    SynthRain$K = PETresults$K;
    PETProb = PETresults$PETProb;
    rainProb = rainSeriesResults$DaysProb;
    SynthRainEP = SynthRain %>% arrange(year, dayIndex);
    #SynthRainS %>% filter(rain > 0) %>% group_by(year) %>% summarise(s = sum(rain), n = n()) %>% ungroup %>% summarise_all("mean")

    #sedom
    IMSRain = GetImsRain(station = 337000, stationEvap = 337000);
    rainSeriesResults = GenerateSeries(NumOfSeries = 500, IMSRain = IMSRain, AnuualRain = array$annual[i], WetDays = array$days[i])

    PETresults = PETGen(rainSeriesResults$SynthRain, IMSRain, 30);
    SynthRain = rainSeriesResults$SynthRain;
    SynthRain$PET = PETresults$SynthPET;
    SynthRain$K = PETresults$K;
    PETProb = PETresults$PETProb;
    rainProb = rainSeriesResults$DaysProb;
    SynthRainSP = SynthRain %>% arrange(year, dayIndex);

    clusterExport(cl, "SynthRainEP")
    clusterExport(cl, "SynthRainSP")


    #rain sens ---

    results$T1.10 = c(results$T1.10, parLapply(cl, 1:5, fun = function(X) CalcGypsum(SynthRainE, duration = 20000, Depth = 150, distribution = T, rainCa = 12, rainSO4 = 12)))
    results$T1.10 = c(results$T1.10, parLapply(cl, 1, fun = function(X) CalcGypsum(SynthRainEP, duration = 20000, Depth = 150, distribution = T, rainCa = 15, rainSO4 = 15)))
    results$T1.10 = c(results$T1.10, parLapply(cl, 1, fun = function(X) CalcGypsum(SynthRainEP, duration = 20000, Depth = 150, distribution = T, rainCa = 20, rainSO4 = 20)))
    results$T1.10 = c(results$T1.10, parLapply(cl, 1, fun = function(X) CalcGypsum(SynthRainEP, duration = 20000, Depth = 150, distribution = T, rainCa = 25, rainSO4 = 25)))
    results$T1.10 = c(results$T1.10, parLapply(cl, 1, fun = function(X) CalcGypsum(SynthRainEP, duration = 20000, Depth = 150, distribution = T, rainCa = 30, rainSO4 = 30)))

    results$T1.9 = c(results$T1.9, parLapply(cl, 1:10, fun = function(X) CalcGypsum(SynthRainEP, duration = 11000, Depth = 150, distribution = T)))
    results$zel11 = c(results$zel11, parLapply(cl, 1:10, fun = function(X) CalcGypsum(SynthRainSP, duration = 10300, Depth = 150, distribution = T)))
}

save(results, file = "resultsListCalib.RData")

future::plan(sequential)
resultsTable = bind_rows(
   RMSD.T.10 = RectanglingResults(results$T1.10 , c(T1.10Observed %>% pull(correctedMean))) %>% mutate(profile = "T1.10", isHolocene = T),
    RMSD.T1.9 = RectanglingResults(results$T1.9 , c(T1.9Observed %>% pull(correctedMean))) %>% mutate(profile = "T1.9", isHolocene = T),
    RMSD.Zel11 = RectanglingResults(results$zel11 , c(Zel11Observed %>% pull(correctedMean))) %>% mutate(profile = "zel11", isHolocene = T),
#RMSD.T1.1 = RectanglingResults(results$T1.1, c(T1.1Observed %>% pull(correctedMean))) %>% mutate(profile = "T1.1", isHolocene = F),
#RMSD.T2.1 = RectanglingResults(results$T2.1, c(T2.1Observed %>% pull(correctedMean))) %>% mutate(profile = "T2.1", isHolocene = F),
#RMSD.zel1 = RectanglingResults(results$zel1, c(Zel1Observed %>% pull(correctedMean))) %>% mutate(profile = "zel1", isHolocene = F)
)
future::plan(strategy = sequential)
resultsTable = resultsTable %>% mutate(optimal = ifelse(FC == opt.FC & AETF == opt.AETF & sulfate == opt.sulfate & dustFlux == opt.dust & WP == opt.WP, T, F)) %>%
    mutate(region = ifelse(profile %in% c("zel11", "zel1"), "zeelim", "shehoret"))


## sens test for weather series---
#%>% group_by(scale = round(scale,1), shape = round(shape,1), region)%>% summarise_if(is.numeric,median) 
sensTest = resultsTable %>% filter(dustFlux == 6, FC == opt.FC, AETF == opt.AETF, WP == opt.WP, sulfate == 13) %>% mutate(meanDay = (AnnualRain / rainDays), Gypsum_depth = PeakDepth, Max_Concentration = PeakConc, Total_concentration = total) %>%
    gather("target", "value", Gypsum_depth, Max_Concentration, Total_concentration)
sensTest %>% filter(target == "Gypsum_depth") %>% ggplot(aes(x = AnnualRain, y = value, color = rainDays)) + geom_point() + facet_grid(region ~ target) + scale_color_gradientn(colours = rainbow(5, rev = T)) + labs(x = "Annual rain [mm]", y = "depth [cm]", color = "#  annual\nrain days", size = "#  annual\nrain days") + scale_y_reverse() + geom_point(aes(x = ObsRain, y = ObsDepth), color = "black", size = 10)

sensTest %>% filter(target != "Gypsum_depth") %>% ggplot(aes(x = AnnualRain, y = value, color = (rainDays))) + geom_point() + facet_wrap(region ~ target, scales = "free_y") + scale_color_gradientn(colours = rainbow(5, rev = T)) + labs(x = "Annual rain [mm]", y = "gypsum[meq/100g soil]", color = "#  annual\nrain days", size = "#  annual\nrain days")
sensTest %>% filter(target != "Gypsum_depth") %>% ggplot(aes(x = shape, y = value, color = (AnnualRain))) + geom_point() + facet_wrap(region ~ target, scales = "free_y") + scale_color_gradientn(colours = rainbow(5, rev = T)) + labs()


sensTest %>% filter(target!= "Gypsum_depth")  %>% ggplot(aes(x = exc, y = value, color = NULL)) + geom_point(size = 2) + facet_wrap(target~profile, scales = "free_y") + scale_color_gradientn(colours = rainbow(5, rev = T)) + labs() +coord_cartesian(xlim = c(50,100))
#---


#print gypsum path
results$T1.10[1] %>% map_dfc(~tibble(.x$YearGyp, .x$YearSulfate, .x$YearCa, .x$rainStat$annual)) %>% mutate_all(~MovingAvarage(.x, 100, 100)) %>% rowid_to_column(var = "year") %>% gather("profile", "totalGypsum", - year) %>% ggplot(aes(year, totalGypsum, color = profile)) + coord_cartesian(xlim = c()) + geom_point()
results$T1.10 %>% map_dfc(~tibble(.x$YearGyp)) %>% mutate_all(~MovingAvarage(.x, 100, 100)) %>% rowid_to_column(var = "year") %>% gather("profile", "totalGypsum", - year) %>% ggplot(aes(year, totalGypsum, color = profile))+ geom_point() 

bla = results$T1.10[c(1)] %>% map_dfc(~tibble(.x$rainStat, gyp = .x$YearGyp, maxgyp = .x$YearMaxGyp, so4 = .x$YearSulfate, ca = .x$YearCa)) %>% mutate(gypagg = gyp - lag(gyp), caagg = ca - lag(ca), so4agg = so4 - lag(so4))
bla %>% dplyr::select(gyp,gyp1,year,PET,PET1) %>% mutate_all( ~ MovingAvarage(.x, 100)) %>% gather("profile", "value", - year) %>%
ggplot(aes(year, value, color = profile)) + geom_point() +labs(y = "[mm] / [meq]")

x=bla %>% filter(gypagg>0) %>% pull(gypagg) 
qnt <- quantile(x, probs = c(.25, .75), na.rm = T)
caps <- quantile(x, probs = c(.05, .95), na.rm = T)
H <- 1.5 * IQR(x, na.rm = T)
bla %>% filter(gypagg > (qnt[1] - H), gypagg < (qnt[2] + H)) %>% ggplot(aes(annual, gypagg)) + geom_point() + coord_cartesian(ylim = c()) + scale_color_gradientn(colours = rainbow(5, rev = T)) + labs(y = "annual gypsum percipitation [meq/100 gr soil]", x = "annual rain [mm]")
bla %>% filter(gypagg < 0) %>% ggplot(aes(annual, gypagg)) + geom_point() + coord_cartesian(ylim = c()) + scale_color_gradientn(colours = rainbow(5, rev = T)) + labs(y = "annual gypsum percipitation [meq/100 gr soil]", x = "annual rain [mm]")
bla %>% ggplot(aes(annual , gypagg, color = n)) + geom_point() + coord_cartesian(ylim = c()) + scale_color_gradientn(colours = rainbow(5, rev = T)) 