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

    RainArr = 24;
DustArr = c(1,10,20,30);
repetition = 1:60;
PETseq = seq(1200, 2400, by = 400)

rainDustArray = tibble(crossing(RainArr, DustArr,PETseq, repetition));

    annual = seq(20, 100, by = 10);
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
    wantedDays = array$days[i];
    wantedRain = array$annual[i];
     #rainDustArray = success %>% filter(WRain == array$days[i], Wdays == array$annual[i]) %>% ungroup() %>% dplyr::select(PETseq = Wpet, sulfate,  dustFlux) %>% unique() %>% dplyr::slice(rep(row_number(), 10))
   # if (nrow(rainDustArray) == 0) {next; }
    #eilat
    SynthRainEP = GenerateSeries(station = 347700, stationEvap = 347704, NumOfSeries = 800, AnuualRain = array$annual[i], WetDays = array$days[i])

    #sedom
   # SynthRainSP = GenerateSeries(station = 337000, stationEvap = 337000, NumOfSeries = 500, AnuualRain = array$annual[i], WetDays = array$days[i])

    clusterExport(cl, list("SynthRainEP","wantedDays","wantedRain"))
    #clusterExport(cl, "rainDustArray")
    rm(SynthRainEP)
    #SynthRainEP %>% filter(rain > 0) %>% group_by(year) %>% summarise(s = sum(rain), n = n()) %>% ungroup %>% summarise_all("mean")        

    #rain sens ---
    print("sendind to calculation: 1")
    shehoret = parLapply(cl, 1:nrow(rainDustArray), fun = function(X) CalcGypsum(SynthRainE, SynthRainEP, duration = T2.1$AvgAge[1], Depth = 150, rainSO4 = rainDustArray$RainArr[X], dustFlux = rainDustArray$DustArr[X], PETFactor = getPETfactor(347700, rainDustArray$PETseq[X]), WRain = wantedRain, Wdays = wantedDays, Wpet = rainDustArray$PETseq[X]))
    print("2")
    #results$Qa2 = c(results$Qa2, parLapply(cl, 1:nrow(rainDustArray), fun = function(X) CalcGypsum(SynthRainE, SynthRainEP, duration = T2.10$AvgAge[1], Depth = 150, rainSO4 = rainDustArray[X, 1], dustFlux = rainDustArray[X, 2], AETFactor = rainDustArray[X, 3])))
    # print("2")
    #zeelim = parLapply(cl, 1:nrow(rainDustArray), fun = function(X) CalcGypsum(SynthRainS, SynthRainSP, duration = Zel1$AvgAge[1], Depth = 150, rainSO4 = rainDustArray$RainArr[X], dustFlux = rainDustArray$DustArr[X], PETFactor = getPETfactor(337000, rainDustArray$PETseq[X]), WRain = wantedRain, Wdays = wantedDays, Wpet = rainDustArray$PETseq[X]))

    results = c(results, zeelim, shehoret)
    # results$talus1 = c(results$talus1, parLapply(cl, 1:nrow(rainDustArray), fun = function(X) CalcGypsum(SynthRainS, SynthRainSP, duration = Talus1Observed$AvgAge[1], Depth = 150, rainSO4 = rainDustArray[X, 1], dustFlux = rainDustArray[X, 2], AETFactor = rainDustArray[X, 3])))

    toc()
    cat(paste(i, " out of ", nrow(array), " ", Sys.time(), "\n"), file = "C:\\Users\\liorst\\Google Drive\\outfile.txt", append = T)
    save(results, file = "resultsListCalibRainPleistAfterSdomCorrection4.RData")

}

save(results, resultsTable, file = "resultsListCalibRainAllScenarios.RData")

resultsTable = bind_rows(resultsTable,a)
results = c(results1,results)
resultsTable = RectanglingResults(results, T1.1) %>% mutate(isHolocene = F) %>% mutate(surface = ifelse(duration == 62500, "Qa1", "Terrace1"), region = ifelse(duration == 62500, "Shehoret", "Zeelim"))
resultsTable = bind_cols(resultsTable, crossing(array, 1:2, rainDustArray))

surfaces = tibble(surface = c("Qa1","Terrace1"), top = c(6,30),bottom = c(20,70), Obstotal = c(70,17.3), sd = c(3.3,3),obsDays = c(9.6,15),obsMAR = c(17.6,39),obsPET = c(2100,2300))
surfaces = resultsTable %>% left_join(surfaces, by = "surface")%>% mutate(PET = PET/100,Wpet=Wpet/100)
surfaces = surfaces %>% mutate(inObsRange = (total <= Obstotal * 1.1 & total >= Obstotal * 0.9 & PeakDepth <= bottom & PeakDepth >= top)) %>% group_by(region, WRain, Wdays, Wpet, sulfate, dustFlux) %>% mutate(id = cur_group_id())

scenarioTable = surfaces %>% group_by(region, WRain, Wdays, Wpet, sulfate, dustFlux,duration) %>% summarise(.groups = "keep", id = cur_group_id(), AnnualRain = mean(AnnualRain), rainDays = mean(rainDays), PET = mean(PET), timesSimulated = n(), timesScored = sum(inObsRange), fitRate = timesScored / timesSimulated * 100, obsMAR = mean(Obstotal), obsDays = mean(obsDays), obsPET = mean(obsPET)) %>% filter() %>% group_by(duration) %>% mutate(weight = timesSimulated / sum(timesSimulated), normalScor = weight * fitRate)

ggarrange(
scenarioTable %>% filter(fitRate > 0) %>% group_by(WRain, sulfate, region) %>% summarise_all(mean) %>% ggplot(aes(WRain, sulfate, fill = fitRate)) + geom_raster(interpolate = F) + scale_fill_distiller(palette = "Greys", direction = 1) + facet_wrap(. ~ region)+theme(legend.position = "none"),
scenarioTable %>% filter(fitRate > 0) %>% group_by(WRain, Wdays, region) %>% summarise_all(mean) %>% ggplot(aes(WRain, Wdays, fill = fitRate)) + geom_raster(interpolate = F) + scale_fill_distiller(palette = "Greys", direction = 1) + facet_wrap(. ~ region) + theme(legend.position = "none"),
scenarioTable %>% filter(fitRate > 0) %>% group_by(WRain, Wpet, region) %>% summarise_all(mean) %>% ggplot(aes(WRain, Wpet, fill = fitRate)) + geom_raster(interpolate = F) + scale_fill_distiller(palette = "Greys", direction = 1) + facet_wrap(. ~ region) + theme(legend.position = "none"),
scenarioTable %>% filter(fitRate > 0) %>% group_by(WRain, dustFlux, region) %>% summarise_all(mean) %>% ggplot(aes(WRain, dustFlux, fill = fitRate)) + geom_raster(interpolate = F) + scale_fill_distiller(palette = "Greys", direction = 1) + facet_wrap(. ~ region) + theme(legend.position = "none"),
scenarioTable %>% filter(fitRate > 0) %>% group_by(sulfate, Wpet, region) %>% summarise_all(mean) %>% ggplot(aes(sulfate, Wpet, fill = fitRate)) + geom_raster(interpolate = F) + scale_fill_distiller(palette = "Greys", direction = 1) + facet_wrap(. ~ region) + theme(legend.position = "none"),
scenarioTable %>% filter(fitRate > 0) %>% group_by(sulfate, dustFlux, region) %>% summarise_all(mean) %>% ggplot(aes(sulfate, dustFlux, fill = fitRate)) + geom_raster(interpolate = F) + scale_fill_distiller(palette = "Greys", direction = 1) + facet_wrap(. ~ region) + theme(legend.position = "none"),
scenarioTable %>% filter(fitRate > 0) %>% group_by(sulfate, Wdays, region) %>% summarise_all(mean) %>% ggplot(aes(sulfate, Wdays, fill = fitRate)) + geom_raster(interpolate = F) + scale_fill_distiller(palette = "Greys", direction = 1) + facet_wrap(. ~ region) + theme(legend.position = "none"),
scenarioTable %>% filter(fitRate > 0) %>% group_by(Wdays, dustFlux, region) %>% summarise_all(mean) %>% ggplot(aes(Wdays, dustFlux, fill = fitRate)) + geom_raster(interpolate = F) + scale_fill_distiller(palette = "Greys", direction = 1) + facet_wrap(. ~ region) + theme(legend.position = "none"),
scenarioTable %>% filter(fitRate > 0) %>% group_by(Wdays, Wpet, region) %>% summarise_all(mean) %>% ggplot(aes(Wdays, Wpet, fill = fitRate)) + geom_raster(interpolate = F) + scale_fill_distiller(palette = "Greys", direction = 1) + facet_wrap(. ~ region) + theme(legend.position = "none"),
scenarioTable %>% filter(fitRate > 0) %>% group_by(Wpet, dustFlux, region) %>% summarise_all(mean) %>% ggplot(aes(Wpet, dustFlux, fill = fitRate)) + geom_raster(interpolate = F) + scale_fill_distiller(palette = "Greys", direction = 1) + facet_wrap(. ~ region) + theme(legend.position = "none")
, ncol = 5)

scenarioTable %>% filter((region == "Shehoret" & fitRate > 70 & WRain == 40 & sulfate <= 22 & sulfate >= 22) | (region == "Zeelim" & fitRate > 50 & WRain == 60 & sulfate == 10)) %>% group_by(sulfate,WRain,Wdays, dustFlux, region) %>% summarise(n = n(), fitRate = mean(fitRate), .groups = "keep") %>% ggplot(aes(dustFlux, Wdays, fill = fitRate)) + geom_raster(interpolate = F) + scale_fill_distiller(palette = "Blues", direction = 1) + facet_wrap(. ~ paste(region,": ", " sulfate: :",sulfate," Annual Rain: ", WRain)) + theme() + coord_cartesian()
scenarioTable %>% filter(region == "Shehoret", fitRate > 75, WRain == 40) %>% ggplot(aes(sulfate, rainDays, color = fitRate)) + geom_point(size = 4) + scale_color_distiller(palette = "Blues", direction = 1) + facet_wrap(. ~ region) + theme()
bla = surfaces %>% group_by(duration, annual, days, PETseq, sulfate, dustFlux) %>% summarise(timesSimulated = n(), timesScored = sum(inObsRange), fitRate = timesScored / timesSimulated * 100)

bla = surfaces %>% dplyr::select(c(1:5, 9, 17, 30:32)) %>% group_by(duration, raingroup, dayGroup, PETGroup, sulfate, dustFlux) %>% summarise(n = n(), group_no = paste(duration, raingroup, dayGroup, PETGroup, sulfate, dustFlux), row_number()) %>% filter(n < 40)

success = scenarioTable %>% filter(fitRate > 50)
#plot map with all points for initial range
E = surfaces %>% filter() %>% group_by(total = round(total, 2), PeakDepth = round(PeakDepth, 2), region, surface) %>% summarise_all(mean) %>% ggplot(aes(total, PeakDepth, color = dustFlux, ymin = bottom, ymax = top, xmin = Obstotal * 0.9, xmax = Obstotal * 1.1)) + geom_point() + facet_wrap(. ~ paste(region, ": ", surface, " ", duration, " years"), scales = "free") + geom_rect(color = "black", fill = NA, size = 2) + scale_y_reverse() + labs(x = "Total gypsum conc. [meq/100 gr soil]", y = "Gypsic horizon depth [cm]") + scale_color_distiller(palette = "PuBuGn", direction = 1)
#plot box plot showing how sulfate and rain are the most sensitive
boxplotstable = scenarioTable %>% dplyr::select(sulfate, dustFlux, Wdays, WRain, region, Wpet, fitRate) %>% gather("key", "value", - region, - fitRate, - duration)
ggplot(boxplotstable %>% filter(fitRate > 85), aes(y = value, x = key, color = region)) + geom_boxplot(outlier.color = NA,position = position_dodge(0.9)) + geom_point(position = position_dodge(0.9)) + scale_y_continuous(breaks = scales::extended_breaks(20)) + theme(axis.title.x = element_blank()) + scale_fill_brewer(type = "qual", palette = 2) + scale_color_brewer(type = "qual", palette = 2) + stat_summary(fun.y = min, fun.ymax = length, geom = "text", aes(label = ..ymax..), vjust = 2, position = position_dodge(0.9)) + coord_cartesian(ylim = c()) +
    stat_summary(data = boxplotstable, fun.min = min,fun = min,fun.max = max, geom = "crossbar", position = position_dodge(0.9),fatten = 0,linetype = "blank", fill = "gray",alpha = 0.5)


scenarioTable %>% dplyr::select(sulfate, dustFlux, Wdays, WRain, region, Wpet, fitRate) %>% gather("key", "value", - region, - fitRate,duration) %>% ggplot(aes(y = value, x = key, color = region, alpha = factor(fitRate > 80))) + geom_boxplot(outlier.color = NA) + geom_point(position = position_dodge(0.75)) + scale_y_continuous(breaks = scales::extended_breaks(20)) + theme(axis.title.x = element_blank()) + scale_fill_brewer(type = "qual", palette = 2) + scale_color_brewer(type = "qual", palette = 2) + stat_summary(fun = median, fun.min = min, fun.max = max, geom = "pointrange", color = "red") + coord_cartesian(ylim = c(0, 150))
a = scenarioTable %>% dplyr::select(sulfate, dustFlux, Wdays, WRain, region, Wpet, fitRate) %>% group_by(region, fitRate = if_else(fitRate > 85, "T", "F")) %>% summarise_all(~length(unique(.))) %>% group_by(region) %>% pivot_longer(cols = c(-region, - fitRate)) %>% pivot_wider(names_from = fitRate) %>% mutate(range = T / F)
 a %>% ggplot(aes(x = name, y = range, fill = region)) + geom_col(position = position_dodge(0.75)) + scale_fill_brewer(type = "qual", palette = 2)
surfaces %>% dplyr::select(surface, duration, sulfate, dustFlux, Wdays, WRain, region, Wpet, inObsRange) %>% gather("key", "value", - surface, - duration, - inObsRange) %>% ggplot(aes(y = value, x = key, color = surface, alpha = factor(inObsRange))) + scale_alpha_discrete(range = c(1, 0.5))  + geom_point(position = position_dodge(0.75)) + scale_y_continuous(breaks = scales::extended_breaks(20)) + theme(axis.title.x = element_blank()) + scale_fill_brewer(type = "qual", palette = 2) + scale_color_brewer(type = "qual", palette = 2) + stat_summary(fun.y = min, fun.ymax = length, geom = "text", aes(label = ..ymax..), vjust = 2, position = position_dodge(0.75))
surfaces %>% filter(region == "Shehoret") %>% dplyr::select(surface, duration, sulfate, dustFlux, AnnualRain, rainDays, region, PET, inObsRange) %>% gather("key", "value", - surface, - duration, - inObsRange, - region) %>% ggplot(aes(y = value, x = key, color = factor(inObsRange))) + geom_point(size = 4) + scale_alpha_discrete(range = c(0.001,1))

scenarioTable %>% filter(fitRate > 80) %>% mutate(sulfate = sulfate / opt.sulfate, dustFlux = dustFlux / opt.dust, Wdays = Wdays / obsDays, WRain = WRain / obsMAR, Wpet = Wpet / obsPET) %>%
    dplyr::select(sulfate, dustFlux, Wdays, WRain, region, Wpet) %>% gather("key", "value", - region) %>% ggplot(aes(y = value, x = key, color = region)) + geom_boxplot(outlier.colour = NA) + geom_point(position = position_dodge(0.75)) + scale_y_continuous(breaks = scales::extended_breaks(20)) + theme(axis.title.x = element_blank()) + scale_fill_brewer(type = "qual", palette = 2) + scale_color_brewer(type = "qual", palette = 2) + labs(y = "ratio from ref [-]") + stat_summary(fun = min, fun.max = length, geom = "text", aes(label = ..ymax..), vjust = 2, position = position_dodge(0.75)) + coord_cartesian(ylim = c(0, 8)) + geom_hline(aes(yintercept = 1), linetype = 2)


#get range of gysic horizons
observedProfiles %>% filter(AvgAge > 50000 & SiteName != "Sheoheret 4" & SiteName != "Shehoret T1-1") %>% group_by(SiteName, correctedMean) %>% summarise(top = min(top), bottom = max(bottom)) %>% slice(n()) %>% ungroup() %>% summarise(min(top), max(bottom))
observedProfiles %>% filter() %>% group_by(SiteName) %>% summarise(correctedMeanTotal = sumGypsum(correctedMean, 5))
observedProfiles %>% filter(AvgAge == 22900) %>% group_by(SiteName) %>% summarise(correctedMeanTotal = sumGypsum(correctedMean, 5)) %>% pull(correctedMeanTotal) %>% mean

scenarioTable1 %>% filter(fitRate > 0, region == "Zeelim") %>% group_by(WRain, Wdays) %>% mutate(name = paste(WRain, ", ", Wdays)) %>% ggplot(aes(x = name, y = sulfate, color = fitRate)) + geom_point(size = 4) + scale_color_gradient2(low = "white", high = "blue")

ggarrange(a, b, c, d, E,nrow = 1)