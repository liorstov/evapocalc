require(parallel)
load("synth.RData")
cl <- makeCluster(4, type = "PSOCK")

clusterExport(cl, c("Zel11Observed", "T1.9Observed", "T1.10Observed", "T2.1", "T1.1Observed", "opt.AETF", "opt.dust", "opt.FC", "opt.sulfate", "opt.WP", "Zel1", "Zel2"))

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
 RMSD.T.10 = RectanglingResults(results$T1.10 , T1.10) %>% mutate(profile = "T1.10", isHolocene = T),
    RMSD.T1.9 = RectanglingResults(results$T1.9 , T1.9) %>% mutate(profile = "T1.9", isHolocene = T),
    RMSD.Zel11 = RectanglingResults(results$zel11 , Zel11) %>% mutate(profile = "zel11", isHolocene = T))

resultsTable = resultsTable %>% mutate(optimal = ifelse(FC == opt.FC & AETF == opt.AETF & sulfate == opt.sulfate & dustFlux == opt.dust & WP == opt.WP, T, F)) %>%
    mutate(region = ifelse(profile %in% c("zel11", "zel1"), "zeelim", "shehoret"))

bservedProfiles %>% filter(AvgAge < 15000) %>% group_by(SiteName) %>% summarise(correctedMeanTotal = sumGypsum(correctedMean, 5))
ShehoretHoloConc = 2.55
ZelHoloConc = 5.23


resultsTable = resultsTable %>% mutate(targetConc = ifelse(region == "zeelim", ZelHoloConc, ShehoretHoloConc)) %>% group_by(sulfate, dustFlux) %>% summarise(n = n(), RMSD = sqrt(sum((total - targetConc) ^ 2) / n))
resultsTable %>% ggplot(aes(sulfate, dustFlux, fill = log10(RMSD), z = log10(RMSD))) + geom_raster(interpolate = T) + scale_fill_gradientn(colors = rainbow(4)) + geom_contour() + labs(x = "Sulfate Concentration in rain water [mg/l]", y = "Annual dust flux [mg/m2/yr]", fill = "RMSD\n[meq/100 gr soil]")

save(results, resultsTable, file = "resultsListCalib.RData")
#rectangle = resultsTable %>% filter(profile == "T1.10") %>% GetResultsInRectangle(T1.10Observed, 30)
#rectangle %>% ggplot(aes(total, PeakDepth, color = inObsRange)) + geom_point()
#rectangle %>% filter(inObsRange) %>% dplyr::select(sulfate, dustFlux, AETF, rainDays, AnnualRain, region) %>% gather("key", "value", - region) %>% ggplot(aes(y = value, x = key, color = region)) + geom_boxplot() + geom_point(position = position_dodge(0.75)) + scale_y_continuous(breaks = scales::extended_breaks(20)) + theme(axis.title.x = element_blank()) + scale_fill_brewer(type = "qual", palette = 2) + scale_color_brewer(type = "qual", palette = 2)


#pareto = bind_rows(paretoSh, paretozel)
#paretozel %>% filter(.level < 10) %>% dplyr::select(sulfate, dustFlux, AETF, rainDays, AnnualRain, region) %>% gather("key", "value", - region) %>% ggplot(aes(y = value, x = key, color = region)) + geom_boxplot() + geom_point(position = position_dodge(0.75)) + scale_y_continuous(breaks = scales::extended_breaks(20)) + theme(axis.title.x = element_blank()) + scale_fill_brewer(type = "qual", palette = 2) + scale_color_brewer(type = "qual", palette = 2)
  #%>% dplyr::select(sulfate, dustFlux, AETF, rainDays, AnnualRain, region) %>% gather("key", "value", - region) %>% ggplot(aes(y = value, x = key, color = region)) + geom_boxplot() + geom_point(position = position_dodge(0.75)) + scale_y_continuous(breaks = scales::extended_breaks(20)) + theme(axis.title.x = element_blank()) + scale_fill_brewer(type = "qual", palette = 2) + scale_color_brewer(type = "qual", palette = 2)
#pareto %>% filter(pareto) %>% dplyr::arrange(distanceOrigin)
#paretoSh = resultsTable %>% group_by(sulfate, dustFlux) %>% summarise_if(is.numeric, "mean") %>% ungroup() %>% psel(low(biasDepth) * low(biasConc), top = 1000) %>% mutate(distanceOrigin = abs(biasConc - RMSD))
#paretoSh %>% ggplot(aes(RMSD, biasConc, shape = .level == 1)) + labs(x = "RMSD [meq/100g soil]", y = "Bias [meq/100g soil]", size = "dust flux mg/m2/yr") + geom_point(size = 3) + scale_color_gradientn(colours = rainbow(5, rev = T)) + scale_shape_manual(values = c(1, 19)) + coord_cartesian(xlim = c(0,20))

RectanglingResults = function(res, obs) {

    ObsDepth = GetGypsicHorizonDepth(obs)
    ObsMaxGyp = GetGypsicHorizonConc(obs)

    obs = obs %>% pull(correctedMean)
    ObsGyp = sumGypsum(obs, 5);
    #calculate target function
    rmsdTable = res %>% map_dfr(.f = ~(tibble(duration = .x$duration, AnnualRain = .x$AnnualRain, rainDays = .x$rainDays, sulfate = .x$RSO4, WP = .x$WP, FC = .x$FC, AETF = .x$AETF, dustFlux = .x$DF, PeakConc = max(.x$gypsum), PeakDepth = comp2Depth(which.max(.x$gypsum), .x$thick),
        total = sumGypsum(.x$gypsum, .x$thick), SWDp80 = .x$SWDp80, obsTop = ObsDepth$top, obsBottom = ObsDepth$bottom, ObsGyp, RMSD = rmsd(.x$gypsum, obs), biasConc = biasConc(.x$gypsum, obs, .x$thick)))) %>% rowid_to_column()
    return(rmsdTable)

}