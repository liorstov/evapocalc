require(parallel)
load("synth.RData")
cl <- makeCluster(1, type = "PSOCK")

clusterExport(cl, c("Zel11Observed", "T1.9Observed", "T1.10Observed", "T2.1Observed","T1.1Observed"))

clusterEvalQ(cl, {
    require(tidyverse)
    require(Rcpp)
    require(ggplot2)
    require(reshape2)
    require(zoo)
    require(tictoc)
    Rcpp::sourceCpp('C:/Users/liorst/source/repos/evapocalc/Calcyp/CSM.cpp', verbose = TRUE, rebuild = 0);
    b <<- new(CSMCLASS);
    source("C:/Users/liorst/source/repos/evapocalc/evapostats/Functions.R", encoding = "Windows-1252")
})

clusterExport(cl, "SynthRainE")
clusterExport(cl, "SynthRainS")
tic()

days = 6:20;
annual =  seq(20, 80, by = 10);
array = tibble(crossing(days, annual),mean = annual/days) %>% filter(mean>=3 & mean <= 11)

for (i in 1:nrow(array)) {
    print(array[i,])

    IMSRain = GetImsRain(station = 347700, stationEvap = 347704);
    rainSeriesResults = GenerateSeries(NumOfSeries = 900, IMSRain = IMSRain, AnuualRain = array$annual[i], WetDays = array$days[i]) 

    PETresults = PETGen(rainSeriesResults$SynthRain, IMSRain, 30);
    SynthRain = rainSeriesResults$SynthRain;
    SynthRain$PET = PETresults$SynthPET;
    SynthRain$K = PETresults$K;
    PETProb = PETresults$PETProb;
    rainProb = rainSeriesResults$DaysProb;
    SynthRain = SynthRain %>% arrange(year, dayIndex);

    clusterExport(cl, "SynthRain")
    
    clusterEvalQ(cl, ls())
  
   # results$T2.1 = c(results$T2.1, parLapply(cl, 1:10, fun = function(X) CalcGypsum(SynthRainE,SynthRain, duration = 62100, plotRes = 0, Depth = tail(T2.1Observed$bottom, 1))))
    #results$T1.1 = c(results$T1.1, parLapply(cl, 1:10, fun = function(X) CalcGypsum(SynthRainE,SynthRain, duration = 62400, plotRes = 0, Depth = tail(T1.1Observed$bottom, 1))))
    results$T1.9 = c(results$T1.9, parLapply(cl, 1:10, fun = function(X) CalcGypsum(SynthRainE, duration = 11000, plotRes = 0, Depth = tail(T1.9Observed$bottom, 1))))
    results$T1.10 = c(results$T1.10, parLapply(cl, 1:10, fun = function(X) CalcGypsum(SynthRainE, duration = 13400, plotRes = 0, Depth = tail(T1.10Observed$bottom, 1))))

   results$zel11 = c(results$zel11,parLapply(cl, 1:10, fun = function(X) CalcGypsum(SynthRainS, duration = 10300, plotRes = 0, Depth = tail(Zel11Observed$bottom, 1))))
   
}
toc()
save(results, file = "resultsList.RData")

stopCluster(cl)

future::plan(multiprocess)
resultsTable = bind_rows(
    RMSD.T.10 = RectanglingResults(results$T1.10 , c(T1.10Observed %>% pull(correctedMean))) %>% mutate(profile = "T1.10", isHolocene = T),
    RMSD.T1.9 = RectanglingResults(results$T1.9 , c(T1.9Observed %>% pull(correctedMean))) %>% mutate(profile = "T1.9", isHolocene = T),
    RMSD.Zel11 = RectanglingResults(results$zel11 , c(Zel11Observed %>% pull(correctedMean))) %>% mutate(profile = "zel11", isHolocene = T),
    #RMSD.T1.1 = RectanglingResults(results$T1.1, c(T1.1Observed %>% pull(correctedMean))) %>% mutate(profile = "T1.1", isHolocene = F),
    #RMSD.T2.1 = RectanglingResults(results$T2.1, c(T2.1Observed %>% pull(correctedMean))) %>% mutate(profile = "T2.1", isHolocene = F)
)
future::plan(strategy = sequential)

#---

plotSoilResultsAgg(results$T1.1[[51]], c(T1.1Observed %>% pull(mean)))