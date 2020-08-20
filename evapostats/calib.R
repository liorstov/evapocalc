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

    RainArr = seq(1, 20, by = 1);
    DustArr = c(seq(1, 10, by = 1));
    repetition = 1
    rainDustArray = as.matrix(crossing(RainArr, DustArr, repetition))
})
results = list()
results$T1.10 = list();
results$T1.9 = list();
results$zel11 = list();


clusterExport(cl, "SynthRainE")
clusterExport(cl, "SynthRainS")
tic()

#calibration ---
results$T1.10 = c(results$T1.10, parLapply(cl, 1:nrow(rainDustArray), fun = function(X) CalcGypsum(SynthRainE, duration = 13400, Depth = 150, rainSO4 = rainDustArray[X, 1], dustFlux = rainDustArray[X, 2])))
results$T1.9 = c(results$T1.9, parLapply(cl, 1:nrow(rainDustArray), fun = function(X) CalcGypsum(SynthRainE, duration = 11000, Depth = 150, rainSO4 = rainDustArray[X, 1], dustFlux = rainDustArray[X, 2])))
results$zel11 = c(results$zel11, parLapply(cl, 1:nrow(rainDustArray), fun = function(X) CalcGypsum(SynthRainS, duration = 10300, Depth = 150, rainSO4 = rainDustArray[X, 1], dustFlux = rainDustArray[X, 2])))

resultsTable = bind_rows(
   RMSD.T.10 = RectanglingResults(results$T1.10, c(T1.10Observed %>% pull(correctedMean))) %>% mutate(profile = "T1.10", isHolocene = T),
    RMSD.T1.9 = RectanglingResults(results$T1.9, c(T1.9Observed %>% pull(correctedMean))) %>% mutate(profile = "T1.9", isHolocene = T),
    RMSD.Zel11 = RectanglingResults(results$zel11, c(Zel11Observed %>% pull(correctedMean))) %>% mutate(profile = "zel11", isHolocene = T),
#RMSD.T1.1 = RectanglingResults(results$T1.1, c(T1.1Observed %>% pull(correctedMean))) %>% mutate(profile = "T1.1", isHolocene = F),
#RMSD.T2.1 = RectanglingResults(results$T2.1, c(T2.1Observed %>% pull(correctedMean))) %>% mutate(profile = "T2.1", isHolocene = F),
#RMSD.zel1 = RectanglingResults(results$zel1, c(Zel1Observed %>% pull(correctedMean))) %>% mutate(profile = "zel1", isHolocene = F)
)
resultsTable = resultsTable %>% mutate(optimal = ifelse(FC == opt.FC & AETF == opt.AETF & sulfate == opt.sulfate & dustFlux == opt.dust & WP == opt.WP, T, F)) %>%
    mutate(region = ifelse(profile %in% c("zel11", "zel1"), "zeelim", "shehoret"))

test = resultsTable %>% filter() %>% groupByParam %>% calculatePareto() #%>% filter(!pareto) %>% calculatePareto %>% filter(!pareto) %>% calculatePareto
test %>% filter(pareto) %>% dplyr::select(sulfate, dustFlux) %>% gather() %>% ggplot(aes(y = value, x = key)) + geom_boxplot() + geom_point() + scale_y_continuous(breaks = scales::extended_breaks(20)) + theme(axis.title.x = element_blank())
test %>% filter(pareto) %>% dplyr::select(sulfate, dustFlux) %>% ungroup %>% summarise_all(median)

save(results, file = "resultsListcalib.RData")
