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

 
})
results = list()
results$T1.10 = list();
results$T1.9 = list();
results$zel11 = list();


clusterExport(cl, "SynthRainE")
clusterExport(cl, "SynthRainS")


#sensitivity ---
results$T1.10 = c(results$T1.10, parLapply(cl, 1:length(seq.dustSeq), function(X) CalcGypsumDouble(SynthRainE, duration = 13400, dustFlux = seq.dustSeq[X])))
results$T1.9 = c(results$T1.9, parLapply(cl, 1:length(seq.dustSeq), function(X) CalcGypsumDouble(SynthRainE, duration = 11000, dustFlux = seq.dustSeq[X])))
results$zel11 = c(results$zel11, parLapply(cl, 1:length(seq.dustSeq), function(X) CalcGypsumDouble(SynthRainS, duration = 10300, dustFlux = seq.dustSeq[X])))

results$T1.10 = c(results$T1.10, parLapply(cl, 1:length(seq.FC), function(X) CalcGypsumDouble(SynthRainE, duration = 13400, FieldCapacity = seq.FC[X])))
results$T1.9 = c(results$T1.9, parLapply(cl, 1:length(seq.FC), function(X) CalcGypsumDouble(SynthRainE, duration = 11000, FieldCapacity = seq.FC[X])))
results$zel11 = c(results$zel11, parLapply(cl, 1:length(seq.FC), function(X) CalcGypsumDouble(SynthRainS, duration = 10300, FieldCapacity = seq.FC[X])))

results$T1.10 = c(results$T1.10, parLapply(cl, 1:length(seq.AETF), function(X) CalcGypsumDouble(SynthRainE, duration = 13400, AETFactor = seq.AETF[X])))
results$T1.9 = c(results$T1.9, parLapply(cl, 1:length(seq.AETF), function(X) CalcGypsumDouble(SynthRainE, duration = 11000, AETFactor = seq.AETF[X])))
results$zel11 = c(results$zel11, parLapply(cl, 1:length(seq.AETF), function(X) CalcGypsumDouble(SynthRainS, duration = 10300, AETFactor = seq.AETF[X])))

results$T1.10 = c(results$T1.10, parLapply(cl, 1:length(seq.WP), function(X) CalcGypsumDouble(SynthRainE, duration = 13400, wieltingPoint = seq.WP[X])))
results$T1.9 = c(results$T1.9, parLapply(cl, 1:length(seq.WP), function(X) CalcGypsumDouble(SynthRainE, duration = 11000, wieltingPoint = seq.WP[X])))
results$zel11 = c(results$zel11, parLapply(cl, 1:length(seq.WP), function(X) CalcGypsumDouble(SynthRainS, duration = 10300, wieltingPoint = seq.WP[X])))

results$T1.10 = c(results$T1.10, parLapply(cl, 1:length(seq.rainSeq), function(X) CalcGypsumDouble(SynthRainE, duration = 13400, rainSO4 = seq.rainSeq[X])))
results$T1.9 = c(results$T1.9, parLapply(cl, 1:length(seq.rainSeq), function(X) CalcGypsumDouble(SynthRainE, duration = 11000, rainSO4 = seq.rainSeq[X])))
results$zel11 = c(results$zel11, parLapply(cl, 1:length(seq.rainSeq), function(X) CalcGypsumDouble(SynthRainS, duration = 10300, rainSO4 = seq.rainSeq[X])))

results$T1.10 = c(results$T1.10, parLapply(cl, 1:60, function(X) CalcGypsumDouble(SynthRainE, duration = 13400, wieltingPoint = opt.WP * 0.8, FieldCapacity = opt.FC * 0.8)))
results$T1.9 = c(results$T1.9, parLapply(cl, 1:60, function(X) CalcGypsumDouble(SynthRainE, duration = 11000, wieltingPoint = opt.WP * 0.8, FieldCapacity = opt.FC * 0.8)))
results$zel11 = c(results$zel11, parLapply(cl, 1:60, function(X) CalcGypsumDouble(SynthRainS, duration = 10300, wieltingPoint = opt.WP * 0.8, FieldCapacity = opt.FC * 0.8)))

results$T1.10 = c(results$T1.10, parLapply(cl, 1:60, function(X) CalcGypsumDouble(SynthRainE, duration = 13400, wieltingPoint = opt.WP * 1.2, FieldCapacity = opt.FC * 1.2)))
results$T1.9 = c(results$T1.9, parLapply(cl, 1:60, function(X) CalcGypsumDouble(SynthRainE, duration = 11000, wieltingPoint = opt.WP * 1.2, FieldCapacity = opt.FC * 1.2)))
results$zel11 = c(results$zel11, parLapply(cl, 1:60, function(X) CalcGypsumDouble(SynthRainS, duration = 10300, wieltingPoint = opt.WP * 1.2, FieldCapacity = opt.FC * 1.2)))
#----


save(results, file = "resultsListSoileSens.RData")


stopCluster(cl)

resultsTable = bind_rows(
   RMSD.T.10 = RectanglingResults(results$T1.10 %>% flatten, c(T1.10Observed %>% pull(correctedMean))) %>% mutate(profile = "T1.10", isHolocene = T),
    RMSD.T1.9 = RectanglingResults(results$T1.9 %>% flatten, c(T1.9Observed %>% pull(correctedMean))) %>% mutate(profile = "T1.9", isHolocene = T),
    RMSD.Zel11 = RectanglingResults(results$zel11 %>% flatten, c(Zel11Observed %>% pull(correctedMean))) %>% mutate(profile = "zel11", isHolocene = T),
#RMSD.T1.1 = RectanglingResults(results$T1.1, c(T1.1Observed %>% pull(correctedMean))) %>% mutate(profile = "T1.1", isHolocene = F),
#RMSD.T2.1 = RectanglingResults(results$T2.1, c(T2.1 %>% pull(correctedMean))) %>% mutate(profile = "T2.1", isHolocene = F),
#RMSD.zel1 = RectanglingResults(results$zel1, c(Zel1 %>% pull(correctedMean))) %>% mutate(profile = "zel1", isHolocene = F)
)
resultsTable = resultsTable %>% mutate(optimal = ifelse(FC == opt.FC & AETF == opt.AETF & sulfate == opt.sulfate & dustFlux == opt.dust & WP == opt.WP, T, F)) %>%
    mutate(region = ifelse(profile %in% c("zel11", "zel1"), "zeelim", "shehoret"))

#get values of the reference computation
resultsTable$refPeakConc = resultsTable$PeakConc %>% lead(1)
resultsTable$refPeakDepth = resultsTable$PeakDepth %>% lead(1)
resultsTable$reftotal = resultsTable$total %>% lead(1)

#erase ref computation
senstest = resultsTable %>% filter(rowid %% 2 != 0)

#relative values from ref
senstest = senstest %>% mutate(Gypsum_depth = (PeakDepth - refPeakDepth) / refPeakDepth, Total_concentration = (total - reftotal) / reftotal, WHC = FC - WP);

#build tibble with values for each parameter and calculate the normelise factor
senstest = bind_rows(
   # WHCSens = senstest %>% filter(AETF == opt.AETF, sulfate == opt.sulfate, dustFlux == opt.dust, WP %in% seq.WP) %>% mutate(param = "WHC", change = WHC / opt.WHC, normal = opt.WHC / (WHC - opt.WHC)),
   # WPSens = senstest %>% filter(FC == opt.FC, AETF == opt.AETF, sulfate == opt.sulfate, dustFlux == opt.dust, WP %in% seq.WP) %>% mutate(param = "\u03B8r - residual\nwater content", change = WP / opt.WP, normal = opt.WP / (WP - opt.WP)),
    FCSens = senstest %>% filter(WP == opt.WP, AETF == opt.AETF, sulfate == opt.sulfate, FC %in% seq.FC) %>% mutate(param = "Field Capacity", change = FC / opt.FC, normal = opt.FC / (FC - opt.FC)),
    AETFSens = senstest %>% filter(FC == opt.FC, WP == opt.WP, sulfate == opt.sulfate, dustFlux == opt.dust, AETF %in% seq.AETF) %>% mutate(param = "AET Factor", change = AETF / opt.AETF, normal = opt.AETF / (AETF - opt.AETF)),
    sulfateSens = senstest %>% filter(FC == opt.FC, AETF == opt.AETF, WP == opt.WP, dustFlux == opt.dust, sulfate %in% seq.rainSeq) %>% mutate(param = "Sulfate", change = sulfate / opt.sulfate, normal = opt.sulfate / (sulfate - opt.sulfate)),
    dustFSens = senstest %>% filter(FC == opt.FC, AETF == opt.AETF, sulfate == opt.sulfate, WP == opt.WP, dustFlux %in% seq.dustSeq) %>% mutate(param = "Annual dust flux", change = dustFlux / opt.dust, normal = opt.dust / (dustFlux - opt.dust))
    )

#calculate with normalisation factor
senstest = senstest %>% mutate(normal = abs(normal), Gypsum_depth = Gypsum_depth * normal, Total_concentration = Total_concentration * normal) %>% drop_na()

#senstest$param = factor(senstest$param, levels = c("\u03B8r", "FC", "WHC", "AET.F", "dustFlux", "sulfate"))
#plot
senstest %>% dplyr::select(change,  Gypsum_depth, Total_concentration, Parameter = param,rowid,region) %>% gather("target", "value", - Parameter, - change,-rowid,-region) %>%
    ggplot(aes(y = value, x = factor(change), color = Parameter)) + facet_wrap(target ~ ., nrow = 1) + geom_boxplot(outlier.shape = NA,fill = "gray90") + stat_boxplot(geom = 'errorbar', linetype = 1, width = 0.5, position = position_dodge(0.75)) + geom_hline(aes(yintercept = 0), linetype = "longdash") +
    scale_y_continuous(name = "Rel sensitivity [-]") + theme(strip.text = element_text(size = 30), legend.text = element_text(size = 30), legend.title = element_text(size = 35)) + labs(color = "Parameter", x = "Change from reference")  + scale_x_discrete(labels = function(x) paste0(as.numeric(x) * 100, "%")) + coord_cartesian(ylim = c(-5,10))
