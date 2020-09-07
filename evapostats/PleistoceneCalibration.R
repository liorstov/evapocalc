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
results$T1.1 = list();
results$T2.1 = list();
results$zel1 = list();
results$zel2 = list();



days = seq(5, 30, by = 5);
annual = seq(30, 100, by = 30);
array = tibble(crossing(days, annual), mean = annual / days) %>% filter(mean >= 3 & mean <= 11) %>% arrange(desc(annual))



for (i in 1:nrow(array)) {
    print(paste("iteration: ", i))
    #eilat
    SynthRainEP = GenerateSeries(station = 347700, stationEvap = 347704, NumOfSeries = 800, AnuualRain = array$annual[i], WetDays = array$days[i])

    #sedom
    SynthRainSP = GenerateSeries(station = 337000, stationEvap = 337000, NumOfSeries = 500, AnuualRain = array$annual[i], WetDays = array$days[i])

    #SynthRainSP %>% filter(rain > 0) %>% group_by(year) %>% summarise(s = sum(rain), n = n()) %>% ungroup %>% summarise_all("mean")

    clusterExport(cl, "SynthRainEP")
    clusterExport(cl, "SynthRainSP")
    rm(SynthRainEP, SynthRainSP)

    #rain sens ---

    results$T1.1 = c(results$T1.1, parLapply(cl, 1:nrow(rainDustArray), fun = function(X) CalcGypsum(SynthRainE, SynthRainEP, duration = T1.1Observed$AvgAge[1], Depth = 150, rainSO4 = rainDustArray[X, 1], dustFlux = rainDustArray[X, 2], AETFactor = rainDustArray[X, 3])))
    results$T2.1 = c(results$T2.1, parLapply(cl, 1:nrow(rainDustArray), fun = function(X) CalcGypsum(SynthRainE, SynthRainEP, duration = T2.1Observed$AvgAge[1], Depth = 150, rainSO4 = rainDustArray[X, 1], dustFlux = rainDustArray[X, 2], AETFactor = rainDustArray[X, 3])))
    results$zel1 = c(results$zel1, parLapply(cl, 1:nrow(rainDustArray), fun = function(X) CalcGypsum(SynthRainS, SynthRainSP, duration = Zel1Observed$AvgAge[1], Depth = 150, rainSO4 = rainDustArray[X, 1], dustFlux = rainDustArray[X, 2], AETFactor = rainDustArray[X, 3])))
    results$zel2 = c(results$zel2, parLapply(cl, 1:nrow(rainDustArray), fun = function(X) CalcGypsum(SynthRainS, SynthRainSP, duration = Zel2Observed$AvgAge[1], Depth = 150, rainSO4 = rainDustArray[X, 1], dustFlux = rainDustArray[X, 2], AETFactor = rainDustArray[X, 3])))

    cat(paste(i," out of ",nrow(array)), file = "C:\\Users\\liorst\\Google Drive\\outfile.txt", append = F)
}


resultsTable = bind_rows(
RMSD.T.10 = RectanglingResults(results$T1.1 , c(T1.10Observed %>% pull(correctedMean))) %>% mutate(profile = "T1.10", isHolocene = T),
RMSD.T1.1 = RectanglingResults(results$T2.1, c(T2.1Observed %>% pull(correctedMean))) %>% mutate(profile = "T2.1", isHolocene = F),
RMSD.T2.1 = RectanglingResults(results$zel2, c(T2.1Observed %>% pull(correctedMean))) %>% mutate(profile = "T2.1", isHolocene = F),
RMSD.zel1 = RectanglingResults(results$zel1, c(Zel1Observed %>% pull(correctedMean))) %>% mutate(profile = "zel1", isHolocene = F)
)
resultsTable = resultsTable %>% mutate(optimal = ifelse(FC == opt.FC & AETF == opt.AETF & sulfate == opt.sulfate & dustFlux == opt.dust & WP == opt.WP, T, F)) %>%
    mutate(region = ifelse(profile %in% c("zel11", "zel1"), "zeelim", "shehoret"))

save(results, resultsTable, file = "resultsListCalibRainPleist.RData")

paretoSh = resultsTable %>% filter(profile != "zel1") %>% calculatePareto() 
paretozel = resultsTable %>% filter(profile == "zel1") %>% calculatePareto()
 %>% filter(!pareto) %>% calculatePareto()
%>% calculatePareto %>% filter(!pareto) %>% calculatePareto

rectangle = resultsTable %>% filter(profile == "T2.1") %>% GetResultsInRectangle(T2.1Observed, 10)
rectangle %>% ggplot(aes(total, PeakDepth, color = inObsRange)) + geom_point()
rectangle %>% filter(inObsRange) %>% dplyr::select(sulfate, dustFlux, AETF, rainDays, AnnualRain, region) %>% gather("key", "value", - region) %>% ggplot(aes(y = value, x = key, color = region)) + geom_boxplot() + geom_point(position = position_dodge(0.75)) + scale_y_continuous(breaks = scales::extended_breaks(20)) + theme(axis.title.x = element_blank()) + scale_fill_brewer(type = "qual", palette = 2) + scale_color_brewer(type = "qual", palette = 2)


pareto = bind_rows(paretoSh,paretozel)
paretozel %>% filter(.level<10) %>% dplyr::select(sulfate, dustFlux, AETF, rainDays, AnnualRain, region) %>% gather("key", "value", - region) %>% ggplot(aes(y = value, x = key, color = region)) + geom_boxplot() + geom_point(position = position_dodge(0.75)) + scale_y_continuous(breaks = scales::extended_breaks(20)) + theme(axis.title.x = element_blank()) + scale_fill_brewer(type = "qual", palette = 2) + scale_color_brewer(type = "qual", palette = 2) 
  %>% dplyr::select(sulfate, dustFlux, AETF, rainDays, AnnualRain, region) %>% gather("key", "value", - region) %>% ggplot(aes(y = value, x = key, color = region)) + geom_boxplot() + geom_point(position = position_dodge(0.75)) + scale_y_continuous(breaks = scales::extended_breaks(20)) + theme(axis.title.x = element_blank()) + scale_fill_brewer(type = "qual", palette = 2) + scale_color_brewer(type = "qual", palette = 2)
pareto %>% filter(pareto) %>% dplyr::arrange(distanceOrigin) 
paretoSh = resultsTable %>% filter(profile != "zel1") %>% calculatePareto()
