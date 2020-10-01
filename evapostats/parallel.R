require(parallel)
load("synth.RData")
cl <- makeCluster(4, type = "PSOCK")

clusterExport(cl, c("Zel11Observed", "T1.9Observed", "T1.10Observed", "T2.1", "T1.1Observed", "opt.AETF", "opt.dust", "opt.FC", "opt.sulfate", "opt.WP", "Zel1"))

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

        RainArr = seq(10, 20, by = 1);
    DustArr = c(seq(1, 26, by = 5));
    repetition = 1:5
    rainDustArray = as.matrix(crossing(RainArr, DustArr, repetition))
})
results = list()
results$T1.10 = list();
results$T1.9 = list();
results$zel11 = list();


clusterExport(cl, "SynthRainE")
clusterExport(cl, "SynthRainS")
tic()

days = seq(5, 20, by = 5);
annual =  seq(150, 200, by = 15);
array = tibble(crossing(days, annual),mean = annual/days) %>% filter(mean>=3 & mean <= 11)

for (i in 1:nrow(array)) {
    # eilat
    IMSRain = GetImsRain(station = 347700, stationEvap = 347704);
    rainSeriesResults = GenerateSeries(NumOfSeries = 500, IMSRain = IMSRain, AnuualRain = array$annual[i], WetDays = array$days[i])

    PETresults = PETGen(rainSeriesResults$SynthRain, IMSRain, 30);
    SynthRain = rainSeriesResults$SynthRain;
    SynthRain$PET = PETresults$SynthPET;
    SynthRain$K = PETresults$K;
    PETProb = PETresults$PETProb;
    rainProb = rainSeriesResults$DaysProb;
    SynthRainEP = SynthRain %>% arrange(year, dayIndex);
    #SynthRainS %>% filter(rain > 0) %>% group_by(year) %>% summarise(s = sum(rain), n = n()) %>% ungroup %>% summarise_all("mean")

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
    clusterExport(cl, "SynthRainEP")
    clusterExport(cl, "SynthRainSP")


    #results$T2.1 = c(results$T2.1, parLapply(cl, 1:nrow(rainDustArray), fun = function(X) CalcGypsum(SynthRainE, SynthRainEP, duration = T2.1$AvgAge[1], plotRes = 0, Depth = tail(T2.1$bottom, 1), DustGyp = 0.01, rainSO4 = rainDustArray[X, 1], dustFlux = rainDustArray[X, 2], FieldCapacity = 0.15)))
    #results$zel1 = c(results$zel1, parLapply(cl, 1:nrow(rainDustArray), fun = function(X) CalcGypsum(SynthRainS, SynthRainSP, duration = Zel1$AvgAge[1], plotRes = 0, Depth = tail(Zel1$bottom, 1), DustGyp = 0.01, rainSO4 = rainDustArray[X, 1], dustFlux = rainDustArray[X, 2], FieldCapacity = 0.15)))
    #results$T1.1 = c(results$T1.1, parLapply(cl, 1:10, fun = function(X) CalcGypsum(SynthRainE,SynthRain, duration = 62400, plotRes = 0, Depth = tail(T1.1Observed$bottom, 1))))

    #rain sens ---
    results$T1.10 = c(results$T1.10, parLapply(cl, 1:20, fun = function(X) CalcGypsum(SynthRainEP, duration = 13400, Depth = 200)))
    results$T1.9 = c(results$T1.9, parLapply(cl, 1:20, fun = function(X) CalcGypsum(SynthRainEP, duration = 11000, Depth = 200)))
    results$zel11 = c(results$zel11, parLapply(cl, 1:20, fun = function(X) CalcGypsum(SynthRainSP, duration = 10300, Depth = 200)))
}
    #calibration ---
    #results$T1.10 = c(results$T1.10, parLapply(cl, 1:nrow(rainDustArray), fun = function(X) CalcGypsum(SynthRainE, duration = 13400, Depth = 100, rainSO4 = rainDustArray[X, 1], dustFlux = rainDustArray[X, 2])))
    #results$T1.9 = c(results$T1.9, parLapply(cl, 1:nrow(rainDustArray), fun = function(X) CalcGypsum(SynthRainE, duration = 11000, Depth = 100, rainSO4 = rainDustArray[X, 1], dustFlux = rainDustArray[X, 2])))
    #results$zel11 = c(results$zel11, parLapply(cl, 1:nrow(rainDustArray), fun = function(X) CalcGypsum(SynthRainS, duration = 10300, Depth = 100, rainSO4 = rainDustArray[X, 1], dustFlux = rainDustArray[X, 2])))
    
    #sensitivity ---
    #results$T1.10 = c(results$T1.10, parLapply(cl, 1:length(seq.dustSeq), function(X) CalcGypsumDouble(SynthRainE, duration = 13400, plotRes = 0, Depth = 100, dustFlux = seq.dustSeq[X])))
    #results$T1.9 = c(results$T1.9, parLapply(cl, 1:length(seq.dustSeq), function(X) CalcGypsumDouble(SynthRainE, duration = 11000, plotRes = 0, Depth = 100, dustFlux = seq.dustSeq[X])))
    #results$zel11 = c(results$zel11, parLapply(cl, 1:length(seq.dustSeq), function(X) CalcGypsumDouble(SynthRainS, duration = 10300, plotRes = 0, Depth = 100, dustFlux = seq.dustSeq[X])))

    #results$T1.10 = c(results$T1.10, parLapply(cl, 1:length(seq.FC), function(X) CalcGypsumDouble(SynthRainE, duration = 13400, plotRes = 0, Depth = 100, FieldCapacity = seq.FC[X])))
    #results$T1.9 = c(results$T1.9, parLapply(cl, 1:length(seq.FC), function(X) CalcGypsumDouble(SynthRainE, duration = 11000, plotRes = 0, Depth = 100, FieldCapacity = seq.FC[X])))
    #results$zel11 = c(results$zel11, parLapply(cl, 1:length(seq.FC), function(X) CalcGypsumDouble(SynthRainS, duration = 10300, plotRes = 0, Depth = 100, FieldCapacity = seq.FC[X])))

    #results$T1.10 = c(results$T1.10, parLapply(cl, 1:length(seq.AETF), function(X) CalcGypsumDouble(SynthRainE, duration = 13400, plotRes = 0, Depth = 100, AETFactor = seq.AETF[X])))
    #results$T1.9 = c(results$T1.9, parLapply(cl, 1:length(seq.AETF), function(X) CalcGypsumDouble(SynthRainE, duration = 11000, plotRes = 0, Depth = 100, AETFactor = seq.AETF[X])))
    #results$zel11 = c(results$zel11, parLapply(cl, 1:length(seq.AETF), function(X) CalcGypsumDouble(SynthRainS, duration = 10300, plotRes = 0, Depth = 100, AETFactor = seq.AETF[X])))

    #results$T1.10 = c(results$T1.10, parLapply(cl, 1:length(seq.WP), function(X) CalcGypsumDouble(SynthRainE, duration = 13400, plotRes = 0, Depth = 100, wieltingPoint = seq.WP[X])))
    #results$T1.9 = c(results$T1.9, parLapply(cl, 1:length(seq.WP), function(X) CalcGypsumDouble(SynthRainE, duration = 11000, plotRes = 0, Depth = 100, wieltingPoint = seq.WP[X])))
    #results$zel11 = c(results$zel11, parLapply(cl, 1:length(seq.WP), function(X) CalcGypsumDouble(SynthRainS, duration = 10300, plotRes = 0, Depth = 100, wieltingPoint = seq.WP[X])))

    #results$T1.10 = c(results$T1.10, parLapply(cl, 1:length(seq.rainSeq), function(X) CalcGypsumDouble(SynthRainE, duration = 13400, plotRes = 0, Depth = 100, rainSO4 = seq.rainSeq[X])))
    #results$T1.9 = c(results$T1.9, parLapply(cl, 1:length(seq.rainSeq), function(X) CalcGypsumDouble(SynthRainE, duration = 11000, plotRes = 0, Depth = 100, rainSO4 = seq.rainSeq[X])))
    #results$zel11 = c(results$zel11, parLapply(cl, 1:length(seq.rainSeq), function(X) CalcGypsumDouble(SynthRainS, duration = 10300, plotRes = 0, Depth = 100, rainSO4 = seq.rainSeq[X])))

results$T1.10 = c(results$T1.10, parLapply(cl, 1:30, function(X) CalcGypsumDouble(SynthRainE, duration = 13400, plotRes = 0, Depth = 100, wieltingPoint = opt.WP * 0.8, FieldCapacity=opt.FC*0.8)))
results$T1.9 = c(results$T1.9, parLapply(cl, 1:30, function(X) CalcGypsumDouble(SynthRainE, duration = 11000, plotRes = 0, Depth = 100, wieltingPoint = opt.WP * 0.8, FieldCapacity = opt.FC * 0.8)))
results$zel11 = c(results$zel11, parLapply(cl, 1:30, function(X) CalcGypsumDouble(SynthRainS, duration = 10300, plotRes = 0, Depth = 100, wieltingPoint = opt.WP * 0.8, FieldCapacity = opt.FC * 0.8)))

results$T1.10 = c(results$T1.10, parLapply(cl, 1:30, function(X) CalcGypsumDouble(SynthRainE, duration = 13400, plotRes = 0, Depth = 100, wieltingPoint = opt.WP * 1.2, FieldCapacity = opt.FC * 1.2)))
results$T1.9 = c(results$T1.9, parLapply(cl, 1:30, function(X) CalcGypsumDouble(SynthRainE, duration = 11000, plotRes = 0, Depth = 100, wieltingPoint = opt.WP * 1.2, FieldCapacity = opt.FC * 1.2)))
results$zel11 = c(results$zel11, parLapply(cl, 1:30, function(X) CalcGypsumDouble(SynthRainS, duration = 10300, plotRes = 0, Depth = 100, wieltingPoint = opt.WP * 1.2, FieldCapacity = opt.FC * 1.2)))
#----

   

toc()
save(results, file = "resultsListClimateSens.RData")

stopCluster(cl)

future::plan(multiprocess)
resultsTable = bind_rows(
   RMSD.T.10 = RectanglingResults(results$T1.10 %>% flatten    , c(T1.10Observed %>% pull(correctedMean))) %>% mutate(profile = "T1.10", isHolocene = T),
    RMSD.T1.9 = RectanglingResults(results$T1.9 %>% flatten, c(T1.9Observed %>% pull(correctedMean))) %>% mutate(profile = "T1.9", isHolocene = T),
    RMSD.Zel11 = RectanglingResults(results$zel11 %>% flatten, c(Zel11Observed %>% pull(correctedMean))) %>% mutate(profile = "zel11", isHolocene = T),
    #RMSD.T1.1 = RectanglingResults(results$T1.1, c(T1.1Observed %>% pull(correctedMean))) %>% mutate(profile = "T1.1", isHolocene = F),
    #RMSD.T2.1 = RectanglingResults(results$T2.1, c(T2.1 %>% pull(correctedMean))) %>% mutate(profile = "T2.1", isHolocene = F),
    #RMSD.zel1 = RectanglingResults(results$zel1, c(Zel1 %>% pull(correctedMean))) %>% mutate(profile = "zel1", isHolocene = F)
)
future::plan(strategy = sequential)
resultsTable = resultsTable %>% mutate(optimal = ifelse(FC == opt.FC & AETF == opt.AETF & sulfate == opt.sulfate & dustFlux == opt.dust & WP == opt.WP, T, F)) %>%
    mutate(region = ifelse(profile %in%  c("zel11","zel1"), "zeelim", "shehoret"))
save(resultsTable,results, file = "resultsTableCal.RData")

#---

plotSoilResultsAgg(results$T1.1[[51]], c(T1.1Observed %>% pull(mean)))