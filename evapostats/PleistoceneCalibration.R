require(parallel)
load("synth.RData")
cl <- makeCluster(5, type = "PSOCK")

clusterExport(cl, c("T2.1", "T2.10", "T2.1", "opt.AETF", "opt.dust", "opt.FC", "opt.sulfate", "opt.WP", "Zel1", "Zel2"))

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

    RainArr = seq(10, 35, by = 5);
DustArr = c(seq(1, 31, by = 15));
AetFactorArr = c(seq(1, 1, by = 0.4));
repetition = 1;
rainDustArray = as.matrix(crossing(RainArr, DustArr, AetFactorArr, repetition));

    days = seq(5, 35, by = 5);
annual = seq(30, 150, by = 20);
PET = c(1900,1)
array = tibble(crossing(days, annual,PET), mean = annual / days) %>% filter(mean >= 3 & mean <= 11) 

})

clusterExport(cl, "SynthRainE")
clusterExport(cl, "SynthRainS")
rm(SynthRainE, SynthRainS)
results = list()
results$Terrace1 = list();
results$Qa1 = list();
results$Qa2= list();





for (i in 1:nrow(array)) {
    print(paste("iteration: ", i))
    tic()
    #eilat
    SynthRainEP = GenerateSeries(station = 347700, stationEvap = 347704, NumOfSeries = 800, AnuualRain = array$annual[i], WetDays = array$days[i], PETfactor = array$PET[i])

    #sedom
    SynthRainSP = GenerateSeries(station = 337000, stationEvap = 337000, NumOfSeries = 500, AnuualRain = array$annual[i], WetDays = array$days[i], PETfactor = array$PET[i])

    #SynthRainEP %>% filter(rain > 0) %>% group_by(year) %>% summarise(s = sum(rain), n = n()) %>% ungroup %>% summarise_all("mean")

    clusterExport(cl, "SynthRainEP")
    clusterExport(cl, "SynthRainSP")
    rm(SynthRainEP, SynthRainSP)

    #rain sens ---
    print("sendind to calculation: 1")
    results$Qa1 = c(results$Qa1, parLapply(cl, 1:nrow(rainDustArray), fun = function(X) CalcGypsum(SynthRainE, SynthRainEP, duration = T2.1$AvgAge[1], Depth = 150, rainSO4 = rainDustArray[X, 1], dustFlux = rainDustArray[X, 2], AETFactor = rainDustArray[X, 3])))
    print("2")
    #results$Qa2 = c(results$Qa2, parLapply(cl, 1:nrow(rainDustArray), fun = function(X) CalcGypsum(SynthRainE, SynthRainEP, duration = T2.10$AvgAge[1], Depth = 150, rainSO4 = rainDustArray[X, 1], dustFlux = rainDustArray[X, 2], AETFactor = rainDustArray[X, 3])))
   # print("2")
    results$Terrace1 = c(results$Terrace1, parLapply(cl, 1:nrow(rainDustArray), fun = function(X) CalcGypsum(SynthRainS, SynthRainSP, duration = Zel1$AvgAge[1], Depth = 150, rainSO4 = rainDustArray[X, 1], dustFlux = rainDustArray[X, 2], AETFactor = rainDustArray[X, 3])))
   # results$talus1 = c(results$talus1, parLapply(cl, 1:nrow(rainDustArray), fun = function(X) CalcGypsum(SynthRainS, SynthRainSP, duration = Talus1Observed$AvgAge[1], Depth = 150, rainSO4 = rainDustArray[X, 1], dustFlux = rainDustArray[X, 2], AETFactor = rainDustArray[X, 3])))

    toc()
    cat(paste(i, " out of ", nrow(array), " ", Sys.time(),"\n"), file = "C:\\Users\\liorst\\Google Drive\\outfile.txt",append = T)
    save(results, file = "resultsListCalibRainPleist5.RData")

}



resultsTable = bind_rows(
 RectanglingResults(results$Qa1, T2.1) %>% mutate(surface = "Qa1", isHolocene = F),
 #RectanglingResults(results$Qa2, T1.1) %>% mutate(surface = "Qa2", isHolocene = F),
 RectanglingResults(results$Terrace1, T1.1) %>% mutate(surface = "Terrace1", isHolocene = F)
) %>% mutate(region = ifelse(str_detect("Qa"), "Shehoret","Zeelim"))

surfaces = tibble(surface = c("Qa1","Qa2","Terrace1"), top = c(20,5,20),bottom = c(40,40,35), Obstotal = c(62,6.2,16.6), sd = c(3.3,3,3),obsDays = c(9,9,16),obsMAR = c(17.6,17.6,39),obsPET = c(2100,2100,2300))
surfaces = resultsTable %>% left_join(surfaces, by = "surface")
surfaces = surfaces %>% mutate(inObsRange = total <= Obstotal + sd & total >= Obstotal - sd & PeakDepth <= bottom & PeakDepth >= top) 

surfaces %>% filter() %>% ggplot(aes(total, PeakDepth, color = AnnualRain, ymin = bottom, ymax = top, xmin = Obstotal - sd, xmax = Obstotal + sd)) + geom_point() + facet_wrap(. ~ paste(region,": ",surface, " ", duration, " years"),scales = "free") + geom_rect(color = "black", fill = NA, size = 2) + scale_y_reverse() + labs(x = "Total gypsum conc. [meq/100 gr soil]", y = "Gypsic horizon depth [cm]") + scale_color_gradientn(colors = rainbow(4, rev = 1))
surfaces %>% filter(inObsRange) %>% dplyr::select(sulfate, dustFlux, AETF, rainDays, AnnualRain, surface, duration) %>% gather("key", "value", - surface, - duration) %>% ggplot(aes(y = value, x = key, color = surface)) + geom_boxplot(outlier.color = NA) + scale_y_continuous(breaks = scales::extended_breaks(20)) + theme(axis.title.x = element_blank()) + scale_fill_brewer(type = "qual", palette = 2) + scale_color_brewer(type = "qual", palette = 2)
surfaces %>% filter(inObsRange) %>% mutate(sulfate = sulfate / opt.sulfate, dustFlux = dustFlux / opt.dust, AETF = AETF / opt.AETF, rainDays = rainDays / obsDays, AnnualRain = AnnualRain / obsMAR, PET = PET/obsPET) %>%
    dplyr::select(sulfate, dustFlux, AETF, rainDays, AnnualRain, surface, duration,PET) %>% gather("key", "value", - surface, - duration) %>% ggplot(aes(y = value, x = key, color = surface)) + geom_boxplot(outlier.colour = NA)  + scale_y_continuous(breaks = scales::extended_breaks(20)) + theme(axis.title.x = element_blank()) + scale_fill_brewer(type = "qual", palette = 2) + scale_color_brewer(type = "qual", palette = 2) + labs(y = "diff. from ref [-]")



pareto = bind_rows(paretoSh,paretozel)
paretozel %>% filter(.level<10) %>% dplyr::select(sulfate, dustFlux, AETF, rainDays, AnnualRain, region) %>% gather("key", "value", - region) %>% ggplot(aes(y = value, x = key, color = region)) + geom_boxplot() + geom_point(position = position_dodge(0.75)) + scale_y_continuous(breaks = scales::extended_breaks(20)) + theme(axis.title.x = element_blank()) + scale_fill_brewer(type = "qual", palette = 2) + scale_color_brewer(type = "qual", palette = 2) 
  %>% dplyr::select(sulfate, dustFlux, AETF, rainDays, AnnualRain, region) %>% gather("key", "value", - region) %>% ggplot(aes(y = value, x = key, color = region)) + geom_boxplot() + geom_point(position = position_dodge(0.75)) + scale_y_continuous(breaks = scales::extended_breaks(20)) + theme(axis.title.x = element_blank()) + scale_fill_brewer(type = "qual", palette = 2) + scale_color_brewer(type = "qual", palette = 2)
pareto %>% filter(pareto) %>% dplyr::arrange(distanceOrigin) 
paretoSh = resultsTable %>% filter(profile != "zel1") %>% calculatePareto()

#get range of gysic horizons
observedProfiles %>% filter(AvgAge > 50000 & SiteName != "Sheoheret 4" & SiteName != "Shehoret T1-1") %>% group_by(SiteName, correctedMean) %>% summarise(top = min(top), bottom = max(bottom)) %>% slice(n()) %>% ungroup() %>% summarise(min(top), max(bottom))
observedProfiles %>% filter() %>% group_by(SiteName) %>% summarise(correctedMeanTotal = sumGypsum(correctedMean, 5))  
observedProfiles %>% filter(AvgAge == 22900) %>% group_by(SiteName) %>% summarise(correctedMeanTotal = sumGypsum(correctedMean, 5)) %>% pull(correctedMeanTotal) %>% mean
