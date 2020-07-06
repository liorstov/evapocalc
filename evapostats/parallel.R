require(parallel)
rm(results)
load("synth.RData")
cl <- makeCluster(6, type = "PSOCK")

clusterExport(cl, c("Zel11Observed", "Zel12Observed", "Zel13Observed", "T1.9Observed", "T1.10Observed", "rainDustArray"))

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
rm("SynthRainE")
gc()
clusterEvalQ(cl, ls())
tic()
resultsT1.10.par = parLapply(cl, 1:nrow(rainDustArray), fun = function(X) CalcGypsum(SynthRainE, duration = 13400, plotRes = 0, Depth = tail(T1.10Observed$bottom, 1), dustFlux = rainDustArray[X, 2], rainSO4 = rainDustArray[X, 1]));
resultsT1.9.par = parLapply(cl, 1:nrow(rainDustArray), fun = function(X) CalcGypsum(SynthRainE, duration = 13400, plotRes = 0, Depth = tail(T1.9Observed$bottom, 1), dustFlux = rainDustArray[X, 2], rainSO4 = rainDustArray[X, 1]));
clusterEvalQ(cl, { rm("SynthRainE"); gc(); })
clusterExport(cl, "SynthRainS")
rm("SynthRainS")
gc()
resultsZel11.par = parLapply(cl, 1:nrow(rainDustArray), fun = function(X) CalcGypsum(SynthRainS, duration = 10300, plotRes = 0, Depth = tail(Zel11Observed$bottom, 1), dustFlux = rainDustArray[X, 2], rainSO4 = rainDustArray[X, 1]));
#resultsZel12.par = parLapply(cl, 1:nrow(rainDustArray), fun = function(X) CalcGypsum(SynthRainS, duration = 10300, plotRes = 0, Depth = tail(Zel12Observed$bottom, 1), dustFlux = rainDustArray[X, 2], rainSO4 = rainDustArray[X, 1]));
resultsZel12.par = list();
toc()
#resultsZel13.par = parLapply(cl, 1:nrow(rainDustArray), fun = function(X) CalcGypsum(SynthRainS, duration = 3000, plotRes = 0, Depth = tail(Zel13Observed$bottom, 1), dustFlux = rainDustArray[X, 2], rainSO4 = rainDustArray[X, 1]));
resultsZel13.par = list();
save(list = c("resultsT1.10.par", "resultsT1.9.par", "resultsZel11.par"), file = "test8.RData")
stopCluster(cl)
load(file = "resultsList.RData")
#combine results from different calculations
results = list(s10 = resultsT1.10.par, s9 = resultsT1.9.par, zel11 = resultsZel11.par, zel12 = resultsZel12.par, zel13 = resultsZel13.par) %>% map2(results, append)
save(results, file = "resultsList.RData")

RMSD.T.10 = ArrangeAndCompare(results$s10, c(T1.10Observed %>% pull(mean))) %>% rename(T1.10 = value)
RMSD.T1.9 = ArrangeAndCompare(results$s9, c(T1.9Observed %>% pull(mean))) %>% rename(T1.9 = value)
RMSD.Zel11 = ArrangeAndCompare(results$zel11, c(Zel11Observed %>% pull(mean))) %>% rename(zel11 = value)
#RMSD.Zel12 = ArrangeAndCompare(results$zel12, c(Zel12Observed %>% pull(mean))) %>% rename(zel12 = value)
#RMSD.Zel13 = ArrangeAndCompare(results$zel13, c(Zel13Observed %>% pull(mean))) %>% rename(zel13 = value)
#---
#join all profiles
CV = RMSD.T.10 %>% left_join(RMSD.T1.9, by = c("rain", "dust")) %>% left_join(RMSD.Zel11, by = c("rain", "dust")) %>%
    gather("profile", "value", - rain, - dust) %>% unnest()# %>% rowwise() %>% mutate(minRMSD = joinRMSD(min, comps), meanRMSD = joinRMSD(mean, comps), maxRMSD = joinRMSD(max, comps), normRMSD = joinRMSD(mean, comps) / meanOBS) %>% ungroup() %>% mutate(site = if_else(str_detect(profile, "zel"), "Zel", "SH"))

