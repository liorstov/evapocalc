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
    cppModule <- new(CSMCLASS);
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

    RainArr = seq(8, 22, by = 2);
DustArr = c(1,10,20,30);
repetition = 1:10;
PETseq = seq(1200, 2400, by = 400)

rainDustArray = tibble(crossing(RainArr, DustArr, repetition,PETseq));

    annual = seq(30, 100, by = 10);
days = seq(5, 18, by = 2);
array = tibble(crossing(days, annual), mean = annual / days) %>% filter(mean >= 3 & mean <= 11)

})

clusterExport(cl, "SynthRainE")
clusterExport(cl, "SynthRainS")
rm(SynthRainE, SynthRainS)
results = list()




for (i in 1:nrow(array)) {
    print(paste("iteration: ", i))
    tic()
    #eilat
    SynthRainEP = GenerateSeries(station = 347700, stationEvap = 347704, NumOfSeries = 800, AnuualRain = array$annual[i], WetDays = array$days[i])

    #sedom
    SynthRainSP = GenerateSeries(station = 337000, stationEvap = 337000, NumOfSeries = 500, AnuualRain = array$annual[i], WetDays = array$days[i])

    clusterExport(cl, "SynthRainEP")
    clusterExport(cl, "SynthRainSP")
    rm(SynthRainEP, SynthRainSP)
    #SynthRainEP %>% filter(rain > 0) %>% group_by(year) %>% summarise(s = sum(rain), n = n()) %>% ungroup %>% summarise_all("mean")        

    #rain sens ---
    print("sendind to calculation: 1")
    shehoret = parLapply(cl, 1:nrow(rainDustArray), fun = function(X) CalcGypsum(SynthRainE, SynthRainEP, duration = T2.1$AvgAge[1], Depth = 150, rainSO4 = rainDustArray$RainArr[X], dustFlux = rainDustArray$DustArr[X], PETFactor = getPETfactor(347700, rainDustArray$PETseq[X])))
    print("2")
    #results$Qa2 = c(results$Qa2, parLapply(cl, 1:nrow(rainDustArray), fun = function(X) CalcGypsum(SynthRainE, SynthRainEP, duration = T2.10$AvgAge[1], Depth = 150, rainSO4 = rainDustArray[X, 1], dustFlux = rainDustArray[X, 2], AETFactor = rainDustArray[X, 3])))
    # print("2")
    zeelim = parLapply(cl, 1:nrow(rainDustArray), fun = function(X) CalcGypsum(SynthRainS, SynthRainSP, duration = Zel1$AvgAge[1], Depth = 150, rainSO4 = rainDustArray$RainArr[X], dustFlux = rainDustArray$DustArr[X], PETFactor = getPETfactor(337000, rainDustArray$PETseq[X])))

    results = c(results, zeelim, shehoret)
    # results$talus1 = c(results$talus1, parLapply(cl, 1:nrow(rainDustArray), fun = function(X) CalcGypsum(SynthRainS, SynthRainSP, duration = Talus1Observed$AvgAge[1], Depth = 150, rainSO4 = rainDustArray[X, 1], dustFlux = rainDustArray[X, 2], AETFactor = rainDustArray[X, 3])))

    toc()
    cat(paste(i, " out of ", nrow(array), " ", Sys.time(), "\n"), file = "C:\\Users\\liorst\\Google Drive\\outfile.txt", append = T)
    save(results, file = "resultsListCalibRainPleist1.RData")

}

save(results, resultsTable, file = "resultsListCalibRainPleist.RData")



resultsTable = RectanglingResults(results, T1.1) %>% mutate(isHolocene = F) %>% mutate(surface = ifelse(duration == 62500, "Qa1", "Terrace1"), region = ifelse(duration == 62500, "Shehoret", "Zeelim"))

surfaces = tibble(surface = c("Qa1","Terrace1"), top = c(10,25),bottom = c(40,45), Obstotal = c(63,16.6), sd = c(3.3,3),obsDays = c(9.6,15),obsMAR = c(17.6,39),obsPET = c(2100,2300))
surfaces = resultsTable %>% left_join(surfaces, by = "surface")
surfaces = surfaces %>% mutate(inObsRange = total <= Obstotal * 1.1 & total >= Obstotal * 0.9 & PeakDepth <= bottom & PeakDepth >= top)

bla = surfaces  %>% rowwise() %>% mutate(raingroup = annual[which.min(abs(annual - AnnualRain))],
        dayGroup = days[which.min(abs(days - rainDays))],
        PETGroup = PETseq[which.min(abs(PETseq - PET))])

test = bla %>% group_by(duration, raingroup, dayGroup, PETGroup, sulfate, dustFlux) %>% summarise(timesSimulated = n(), timesScored = sum(inObsRange), fitRate = timesScored / timesSimulated * 100) %>% filter() %>% group_by(duration) %>% mutate(weight = timesSimulated / sum(timesSimulated), normalScor = weight * fitRate)

surfaces %>% filter() %>% ggplot(aes(total, PeakDepth, color = scen, ymin = bottom, ymax = top, xmin = Obstotal * 0.9, xmax = Obstotal * 1.1)) + geom_point() + facet_wrap(. ~ paste(region, ": ", surface, " ", duration, " years"), scales = "free") + geom_rect(color = "black", fill = NA, size = 2) + scale_y_reverse() + labs(x = "Total gypsum conc. [meq/100 gr soil]", y = "Gypsic horizon depth [cm]") + scale_color_gradientn(colors = rainbow(4, rev = 1))
surfaces %>% filter(inObsRange) %>% dplyr::select(sulfate, dustFlux, rainDays, AnnualRain, surface, duration) %>% gather("key", "value", - surface, - duration) %>% ggplot(aes(y = value, x = key, color = surface)) + geom_boxplot(outlier.color = NA) + geom_point(position = position_dodge(0.75)) + scale_y_continuous(breaks = scales::extended_breaks(20)) + theme(axis.title.x = element_blank()) + scale_fill_brewer(type = "qual", palette = 2) + scale_color_brewer(type = "qual", palette = 2) + stat_summary(fun.y = min, fun.ymax = length, geom = "text", aes(label = ..ymax..), vjust = 2, position = position_dodge(0.75)) + coord_cartesian(ylim = c(0, 150))
surfaces %>% filter(inObsRange) %>% dplyr::select(PET, surface, duration) %>% gather("key", "value", - surface, - duration) %>% ggplot(aes(y = value, x = key, color = surface)) + geom_boxplot(outlier.color = NA) + geom_point(position = position_dodge(0.75)) + scale_y_continuous(breaks = scales::extended_breaks(20)) + theme(axis.title.x = element_blank()) + scale_fill_brewer(type = "qual", palette = 2) + scale_color_brewer(type = "qual", palette = 2) + stat_summary(fun.y = min, fun.ymax = length, geom = "text", aes(label = ..ymax..), vjust = 2, position = position_dodge(0.75)) 

surfaces %>% filter(inObsRange) %>% mutate(sulfate = sulfate / opt.sulfate, dustFlux = dustFlux / opt.dust, AETF = AETF / opt.AETF, rainDays = rainDays / obsDays, AnnualRain = AnnualRain / obsMAR, PET = PET / obsPET) %>%
    dplyr::select(sulfate, dustFlux, rainDays, AnnualRain, surface, duration, PET) %>% gather("key", "value", - surface, - duration) %>% ggplot(aes(y = value, x = key, color = surface)) + geom_boxplot(outlier.colour = NA) + geom_point(position = position_dodge(0.75)) + scale_y_continuous(breaks = scales::extended_breaks(20)) + theme(axis.title.x = element_blank()) + scale_fill_brewer(type = "qual", palette = 2) + scale_color_brewer(type = "qual", palette = 2) + labs(y = "ratio from ref [-]") + stat_summary(fun = min, fun.max = length, geom = "text", aes(label = ..ymax..), vjust = 2, position = position_dodge(0.75)) + coord_cartesian(ylim = c(0, 8)) + geom_hline(aes(yintercept = 1), linetype = 2)



#get range of gysic horizons
observedProfiles %>% filter(AvgAge > 50000 & SiteName != "Sheoheret 4" & SiteName != "Shehoret T1-1") %>% group_by(SiteName, correctedMean) %>% summarise(top = min(top), bottom = max(bottom)) %>% slice(n()) %>% ungroup() %>% summarise(min(top), max(bottom))
observedProfiles %>% filter() %>% group_by(SiteName) %>% summarise(correctedMeanTotal = sumGypsum(correctedMean, 5))  
observedProfiles %>% filter(AvgAge == 22900) %>% group_by(SiteName) %>% summarise(correctedMeanTotal = sumGypsum(correctedMean, 5)) %>% pull(correctedMeanTotal) %>% mean
