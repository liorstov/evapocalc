require(parallel)
load("synth.RData")
cl <- makeCluster(4, type = "PSOCK")

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
    opt.sulfate = 10;
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



days = seq(5,30, by = 5);
annual = seq(300,400, by = 10);
array = tibble(crossing(days, annual), mean = annual / days) %>% filter(mean >= 3 & mean <= 11)

for (i in 1:nrow(array)) {
     #eilat
    IMSRain = GetImsRain(station = 347700, stationEvap = 347704);
    rainSeriesResults = GenerateSeries(NumOfSeries = 500, IMSRain = IMSRain, AnuualRain = array$annual[i], WetDays = array$days[i])

    PETresults = PETGen(rainSeriesResults$SynthRain, IMSRain, 30);
    SynthRain = rainSeriesResults$SynthRain;
    SynthRain$PET = PETresults$SynthPET;
    SynthRain$K = PETresults$K;
    PETProb = PETresults$PETProb;
    rainProb = rainSeriesResults$DaysProb;
    SynthRainEP = SynthRain %>% arrange(year, dayIndex);

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
    SynthRainSP %>% filter(rain > 0) %>% group_by(year) %>% summarise(s = sum(rain), n = n()) %>% ungroup %>% summarise_all("mean")

    clusterExport(cl, "SynthRainEP")
    clusterExport(cl, "SynthRainSP")


    #rain sens ---

    results$T1.9 = c(results$T1.9, parLapply(cl, 1:17, fun = function(X) CalcGypsum(SynthRainEP, duration = 11000, Depth = 150, distribution = T)))
    results$zel11 = c(results$zel11, parLapply(cl, 1:17, fun = function(X) CalcGypsum(SynthRainSP, duration = 10300, Depth = 150, distribution = T)))
}

save(results, file = "resultsListCalibRain.RData")

resultsTable = bind_rows(
   RMSD.T.10 = RectanglingResults(results$T1.10 , c(T1.10Observed %>% pull(correctedMean))) %>% mutate(profile = "T1.10", isHolocene = T),
    RMSD.Zel11 = RectanglingResults(results$zel11 , c(Zel11Observed %>% pull(correctedMean))) %>% mutate(profile = "zel11", isHolocene = T),
#RMSD.T1.1 = RectanglingResults(results$T1.1, c(T1.1Observed %>% pull(correctedMean))) %>% mutate(profile = "T1.1", isHolocene = F),
#RMSD.T2.1 = RectanglingResults(results$T2.1, c(T2.1Observed %>% pull(correctedMean))) %>% mutate(profile = "T2.1", isHolocene = F),
#RMSD.zel1 = RectanglingResults(results$zel1, c(Zel1Observed %>% pull(correctedMean))) %>% mutate(profile = "zel1", isHolocene = F)
)
resultsTable = resultsTable %>% mutate(optimal = ifelse(FC == opt.FC & AETF == opt.AETF & sulfate == opt.sulfate & dustFlux == opt.dust & WP == opt.WP, T, F)) %>%
    mutate(region = ifelse(profile %in% c("zel11", "zel1"), "zeelim", "shehoret"))


## sens test for weather series---
#%>% group_by(scale = round(scale,1), shape = round(shape,1), region)%>% summarise_if(is.numeric,median) 
sensTest = resultsTable  %>% mutate(meanDay = (AnnualRain / rainDays), Gypsum_depth = PeakDepth, Max_Concentration = PeakConc, Total_concentration = total) %>%
    gather("target", "value", Gypsum_depth, Max_Concentration, Total_concentration)

sensTest %>% filter(target == "Gypsum_depth") %>% ggplot(aes(x = AnnualRain, y = value, color = rainDays, size = AnnualRain/rainDays)) + geom_point(shape = 1) + facet_wrap(target ~ region, scales = "free_y") + scale_color_gradientn(colours = rainbow(5, rev = T)) + labs() + coord_cartesian(xlim = c()) + scale_y_reverse()
sensTest %>% filter(target != "Gypsum_depth") %>% ggplot(aes(x = AnnualRain, y = value, color = rainDays, size = AnnualRain / rainDays)) + geom_point(shape = 1) + facet_wrap(target ~ region, scales = "free_y") + scale_color_gradientn(colours = rainbow(5, rev = T)) + labs() + coord_cartesian(xlim = c())
#---

