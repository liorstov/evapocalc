require(parallel)
load("synth.RData")
cl <- makeCluster(4, type = "PSOCK")

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
    opt.sulfate = 10;
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



save(results, file = "resultsListSoileSens.RData")

stopCluster(cl)

future::plan(multiprocess)
resultsTable = bind_rows(
   RMSD.T.10 = RectanglingResults(results$T1.10 %>% flatten, c(T1.10Observed %>% pull(correctedMean))) %>% mutate(profile = "T1.10", isHolocene = T),
    RMSD.T1.9 = RectanglingResults(results$T1.9 %>% flatten, c(T1.9Observed %>% pull(correctedMean))) %>% mutate(profile = "T1.9", isHolocene = T),
    RMSD.Zel11 = RectanglingResults(results$zel11 %>% flatten, c(Zel11Observed %>% pull(correctedMean))) %>% mutate(profile = "zel11", isHolocene = T),
#RMSD.T1.1 = RectanglingResults(results$T1.1, c(T1.1Observed %>% pull(correctedMean))) %>% mutate(profile = "T1.1", isHolocene = F),
#RMSD.T2.1 = RectanglingResults(results$T2.1, c(T2.1Observed %>% pull(correctedMean))) %>% mutate(profile = "T2.1", isHolocene = F),
#RMSD.zel1 = RectanglingResults(results$zel1, c(Zel1Observed %>% pull(correctedMean))) %>% mutate(profile = "zel1", isHolocene = F)
)
future::plan(strategy = sequential)
resultsTable = resultsTable %>% mutate(optimal = ifelse(FC == opt.FC & AETF == opt.AETF & sulfate == opt.sulfate & dustFlux == opt.dust & WP == opt.WP, T, F)) %>%
    mutate(region = ifelse(profile %in% c("zel11", "zel1"), "zeelim", "shehoret"))
save(resultsTable, results, file = "resultsTableCal.RData")

#---
senstest = bind_rows(
WPSens = resultsTable %>% filter(FC == opt.FC, AETF == opt.AETF, sulfate == opt.sulfate, dustFlux == opt.dust, WP %in% seq.WP) %>% mutate(param = "\u03B8r", change = WP / opt.WP),
FCSens = resultsTable %>% filter(WP == opt.WP, AETF == opt.AETF, sulfate == opt.sulfate, FC %in% seq.FC) %>% mutate(param = "FC", change = FC / opt.FC),
AETFSens = resultsTable %>% filter(FC == opt.FC, WP == opt.WP, sulfate == opt.sulfate, dustFlux == opt.dust, AETF %in% seq.AETF) %>% mutate(param = "AET.F", change = AETF / opt.AETF),
sulfateSens = resultsTable %>% filter(FC == opt.FC, AETF == opt.AETF, WP == opt.WP, dustFlux == opt.dust, sulfate %in% seq.rainSeq) %>% mutate(param = "sulfate", change = sulfate / opt.sulfate),
dustFSens = resultsTable %>% filter(FC == opt.FC, AETF == opt.AETF, sulfate == opt.sulfate, WP == opt.WP, dustFlux %in% seq.dustSeq) %>% mutate(param = "dustFlux", change = dustFlux / opt.dust)
)
bla = senstest %>% rowid_to_column(var = "id") %>% filter(change == 1) %>% group_by(PeakConc) %>% summarise(mean(FC), n = n(), target = mean(id)) %>% filter(n == 1) %>% pull(target)
senstest = senstest %>% slice(-bla)

senstest = senstest %>% filter(change == 1) %>% dplyr::select(refMaxconc = PeakConc, refDepth = PeakDepth, refTotal = total) %>% summarise_all(median) %>% cbind(senstest) %>% tibble
senstest = senstest %>% mutate(refMaxconc = PeakConc / refMaxconc, refDepth = PeakDepth / refDepth, refTotal = total / refTotal)
save(senstest, file = "senstable.RData")
WP1 = senstest %>% filter(round(change, 1) %in% c(0.8, 1, 1.2)) %>% dplyr::select(change, Max_Concentration = PeakConc, Parameter = param) %>% gather("target", "value", - Parameter, - change) %>%
    ggplot(aes(y = value, x = factor(change), fill = Parameter)) + stat_boxplot(geom = 'errorbar', linetype = 1, width = 0.5, position = position_dodge(0.75)) + geom_boxplot() + geom_point(stat = "summary", fun = "mean", shape = 5, size = 4, position = position_dodge(0.75)) + facet_wrap(target ~ ., scales = "free_y", nrow = 3) +
    scale_y_continuous(name = "Concentration [meq/100g soil]") + theme(strip.text = element_text(size = 30), legend.text = element_text(size = 30), legend.title = element_text(size = 35), axis.title.y = element_text(hjust = 0.9)) + labs(color = "Parameter", x = "change from optimal") + scale_fill_brewer(type = "qual", palette = 7) + scale_x_discrete(labels = function(x) paste0(as.numeric(x) * 100, "%"))
# +geom_boxplot(outlier.shape = NA)
WP2 = senstest %>% filter(round(change, 1) %in% c(0.8, 1, 1.2)) %>% dplyr::select(change, Total_concentration = total, param) %>% gather("target", "value", - param, - change) %>%
    ggplot(aes(y = value, x = factor(change), fill = param)) + stat_boxplot(geom = 'errorbar', linetype = 1, width = 0.5, position = position_dodge(0.75)) + geom_boxplot() + geom_point(stat = "summary", fun = "mean", shape = 5, size = 4, position = position_dodge(0.75)) + facet_wrap(target ~ ., scales = "free_y", nrow = 3) +
    scale_y_continuous(name = "") + theme(strip.text = element_text(size = 30), legend.text = element_text(size = 30), legend.title = element_text(size = 35)) + labs(color = "", x = "change from optimal") + scale_fill_brewer(type = "qual", palette = 7) + scale_x_discrete(labels = function(x) paste0(as.numeric(x) * 100, "%"))
# +geom_boxplot(outlier.shape = NA)
WP3 = senstest %>% filter(round(change, 1) %in% c(0.8, 1, 1.2)) %>% dplyr::select(change, Gypsum_depth = PeakDepth, param) %>% gather("target", "value", - param, - change) %>%
    ggplot(aes(y = value, x = factor(change), fill = param)) + stat_boxplot(geom = 'errorbar', linetype = 1, width = 0.5, position = position_dodge(0.75)) + geom_boxplot() + geom_point(stat = "summary", fun = "mean", shape = 5, size = 4, position = position_dodge(0.75)) + facet_wrap(target ~ ., scales = "free_y", nrow = 3) +
    scale_y_reverse(name = "Depth[cm]", position = "right") + theme(strip.text = element_text(size = 30), axis.text.x = element_text(size = 35), legend.text = element_text(size = 30), legend.title = element_text(size = 35)) + labs(color = " change from \n optimal[%]", x = "Change from optimal") + scale_fill_brewer(type = "qual", palette = 7) + scale_x_discrete(labels = function(x) paste0(as.numeric(x) * 100, "%"))

ggarrange(WP1 + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.key.size = unit(5, "line")),
    WP2 + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position = c("")),
        WP3 + theme(legend.position = c("")))



rainDaysTable = resultsTable %>% filter(profile == "T1.9", FC == opt.FC, AETF == opt.AETF, sulfate == opt.sulfate, dustFlux == opt.dust, WP == opt.WP, !is.na(rainDays)) %>%
    group_by(meanDay = round(AnnualRain / rainDays, 1)) %>% summarise(median = median(total), min = quantile(total, 0.2), max = quantile(total, 0.8)) %>% ggplot(aes(x = meanDay, y = median, ymax = max, ymin = min)) + geom_line() + geom_ribbon(alpha = 0.5)
