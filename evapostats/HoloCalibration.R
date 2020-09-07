require(parallel)
load("synth.RData")
cl <- makeCluster(4, type = "PSOCK")

clusterExport(cl, c("Zel11Observed", "T1.9Observed", "T1.10Observed", "T2.1Observed", "T1.1Observed", "opt.AETF", "opt.dust", "opt.FC", "opt.sulfate", "opt.WP", "Zel1Observed", "Zel2Observed"))

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

    RainArr = seq(10, 30, by = 5);
    DustArr = c(seq(1, 31, by = 15));
    AetFactorArr = c(seq(0.6, 1.4, by = 0.4));
    repetition = 1:1;
    rainDustArray = as.matrix(crossing(RainArr, DustArr, AetFactorArr, repetition));
})

clusterExport(cl, "SynthRainE")
clusterExport(cl, "SynthRainS")
rm(SynthRainE, SynthRainS)
results = list()
results$T1.10 = list();
results$T1.9  = list();
results$zel11 = list();



days = seq(5, 30, by = 5);
annual = seq(30, 100, by = 30);
array = tibble(crossing(days, annual), mean = annual / days) %>% filter(mean >= 3 & mean <= 11) %>% arrange(desc(annual))


#calibration ---
results$T1.10 = c(results$T1.10, parLapply(cl, 1:nrow(rainDustArray), fun = function(X) CalcGypsum(SynthRainE, duration = 13400, Depth = 100, rainSO4 = rainDustArray[X, 1], dustFlux = rainDustArray[X, 2])))
results$T1.9 = c(results$T1.9, parLapply(cl, 1:nrow(rainDustArray), fun = function(X) CalcGypsum(SynthRainE, duration = 11000, Depth = 100, rainSO4 = rainDustArray[X, 1], dustFlux = rainDustArray[X, 2])))
results$zel11 = c(results$zel11, parLapply(cl, 1:nrow(rainDustArray), fun = function(X) CalcGypsum(SynthRainS, duration = 10300, Depth = 100, rainSO4 = rainDustArray[X, 1], dustFlux = rainDustArray[X, 2])))


resultsTable = bind_rows(
 RMSD.T.10 = RectanglingResults(results$T1.10 , T1.10Observed) %>% mutate(profile = "T1.10", isHolocene = T),
    RMSD.T1.9 = RectanglingResults(results$T1.9 , T1.9Observed) %>% mutate(profile = "T1.9", isHolocene = T),
    RMSD.Zel11 = RectanglingResults(results$zel11 , Zel11Observed) %>% mutate(profile = "zel11", isHolocene = T))

resultsTable = resultsTable %>% mutate(optimal = ifelse(FC == opt.FC & AETF == opt.AETF & sulfate == opt.sulfate & dustFlux == opt.dust & WP == opt.WP, T, F)) %>%
    mutate(region = ifelse(profile %in% c("zel11", "zel1"), "zeelim", "shehoret"))


rectangle = resultsTable %>% filter(profile == "T1.10") %>% GetResultsInRectangle(T1.10Observed, 30)
rectangle %>% ggplot(aes(total, PeakDepth, color = inObsRange)) + geom_point()
rectangle %>% filter(inObsRange) %>% dplyr::select(sulfate, dustFlux, AETF, rainDays, AnnualRain, region) %>% gather("key", "value", - region) %>% ggplot(aes(y = value, x = key, color = region)) + geom_boxplot() + geom_point(position = position_dodge(0.75)) + scale_y_continuous(breaks = scales::extended_breaks(20)) + theme(axis.title.x = element_blank()) + scale_fill_brewer(type = "qual", palette = 2) + scale_color_brewer(type = "qual", palette = 2)


pareto = bind_rows(paretoSh, paretozel)
paretozel %>% filter(.level < 10) %>% dplyr::select(sulfate, dustFlux, AETF, rainDays, AnnualRain, region) %>% gather("key", "value", - region) %>% ggplot(aes(y = value, x = key, color = region)) + geom_boxplot() + geom_point(position = position_dodge(0.75)) + scale_y_continuous(breaks = scales::extended_breaks(20)) + theme(axis.title.x = element_blank()) + scale_fill_brewer(type = "qual", palette = 2) + scale_color_brewer(type = "qual", palette = 2)
  %>% dplyr::select(sulfate, dustFlux, AETF, rainDays, AnnualRain, region) %>% gather("key", "value", - region) %>% ggplot(aes(y = value, x = key, color = region)) + geom_boxplot() + geom_point(position = position_dodge(0.75)) + scale_y_continuous(breaks = scales::extended_breaks(20)) + theme(axis.title.x = element_blank()) + scale_fill_brewer(type = "qual", palette = 2) + scale_color_brewer(type = "qual", palette = 2)
pareto %>% filter(pareto) %>% dplyr::arrange(distanceOrigin)
paretoSh = resultsTable %>% filter(profile != "zel1") %>% calculatePareto()
