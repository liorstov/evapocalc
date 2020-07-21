
#load("C:/Users/liorst/Source/Repos/evapocalc/.RData")
load("C:/Users/liorst/Source/Repos/evapocalc/synth.RData")
load(file = "C:/Users/liorst/Source/Repos/evapocalc/resultsList.RData")

require(tidyverse)
require(Rcpp)
require(ggplot2)
require(reshape2)
require(zoo)
require(tictoc)
require(R.matlab)
require(gganimate)
require(pbapply)
require(metR)
require(ggpmisc)
require(akima)
require(rPref)
require(soilptf)
require(furrr)
require(future.apply)

opt.AETF = 1.2;
opt.WP = 0.013;
opt.FC = 0.1;
opt.sulfate = 13;
opt.dust = 6;

theme_set(theme_classic() + theme(legend.key.size = unit(1, "line"), legend.text = element_text(size = 25),
axis.text.x = element_text(size = 28, angle = 43, hjust = 1), title = element_text(size = 20),
strip.text = element_text( size = 20),
axis.text.y = element_text(size = 28),
axis.title.y = element_text(size = 28, margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.text.y.right = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
axis.title.x = element_text(size = 35),
panel.grid.major = element_line(colour = "grey93"),
panel.grid.minor = element_line(colour = "grey93"),
panel.border = element_rect(colour = "black", fill = NA, size = 0))) + guides(colour = guide_legend(override.aes = list(size = 10)));
Sys.setlocale("LC_TIME", "English_Israel.1255");

#set wd 
setwd("C:/Users/liorst/source/repos/evapocalc");


source("evapostats/Functions.R");
source("evapostats/PETGen.R");
source("evapostats/RainGen.R");
loadObsProfiles()

Rcpp::sourceCpp('C:/Users/liorst/source/repos/evapocalc/Calcyp/CSM.cpp', verbose = TRUE, rebuild = 0);
b <<- new(CSMCLASS);

#rain gen ----
#stationElat = 347700;
#stationElatEvap = 347704;
#stationSedom = 337000;
IMSRain = GetImsRain(station = 347700, stationEvap = 347704);
rainSeriesResults = GenerateSeries(NumOfSeries = 900, IMSRain = IMSRain, AnuualRain = 80, WetDays = 25);
PETresults = PETGen(rainSeriesResults$SynthRain, IMSRain,30);
SynthRain = rainSeriesResults$SynthRain;
SynthRain$PET = PETresults$SynthPET;
SynthRain$K = PETresults$K;
PETProb = PETresults$PETProb;
rainProb = rainSeriesResults$DaysProb;
SynthRain = SynthRain %>% arrange(year, dayIndex);

results = plotResults(SynthRain, IMSRain, rainSeriesResults$DaysProb, PETresults$PETProb, 1);

pdf(file = paste("plots/", format(Sys.time(), "%b_%d_%Y_%H%M"), "Results.pdf"), width = 30, height = 16);
print(results);
dev.off()
SynthRain %>% filter(rain > 0) %>% summarise(mean(rain))
SynthRain %>% filter(rain > 0) %>% group_by(year) %>% summarise(n = n()) %>% summarise(mean(n))
#catch observed----

res = tibble(factor = numeric(), wetD = numeric())

resultsT2.1.par = list();
for (days in seq(9, 15)) {
    
    IMSRain = GetImsRain(station = 347700, stationEvap = 347704);
    rainSeriesResults = GenerateSeries(NumOfSeries = 900, IMSRain = IMSRain, AnuualRain = 80, WetDays = days);
    PETresults = PETGen(rainSeriesResults$SynthRain, IMSRain, 30);
    SynthRain = rainSeriesResults$SynthRain;
    SynthRain$PET = PETresults$SynthPET;
    SynthRain$K = PETresults$K;
    PETProb = PETresults$PETProb;
    rainProb = rainSeriesResults$DaysProb;
    SynthRain = SynthRain %>% arrange(year, dayIndex);

    resultsT2.1.par[days] = list()
    resultsT2.1.par[days] = lapply(1:5, FUN = function(X) CalcGypsum(rainHolocene = SynthRainE, rainPleistocene = SynthRain, 62100, plotRes = 0, Depth = tail(T2.1Observed$bottom, 1), rainSO4 = 15, dustFlux = 5, AETFactor = 1.2, random = T, getWD = F));


}

res %>% ggplot(aes(factor, wetD)) + geom_line() +geom_point()

#tests---
result = CalcGypsum(rainHolocene = SynthRainE, rainPleistocene = SynthRain, 62100, plotRes = 0, Depth = tail(T2.1Observed$bottom, 1), rainSO4 = 20, dustFlux = 50, AETFactor = 1.2, random = F, getWD = T);
plotSoilResults(result, T2.1Observed)
bla = result$WD %>% group_by(milen = year %/% 1000) %>% summarise(max = max(meanWD), std = sd(meanWD), mean = mean(meanWD))

ggplot(bla, aes(milen * 1000, mean)) + geom_line(size = 4, color = "blue") +
labs(y = "Mean of seasonal WD[cm]", x = "Soil age [yr]") + scale_y_reverse() + geom_smooth(method = "gam", size = 2,color = "gray2") 

save(resultsT1.1.par, file = "resultsPleis.RData")

kiki1 = SynthRainE %>% filter(year == 1550, dayIndex >= 160 & dayIndex <= 161)


taliRain = readMat("DB\\RainTali.mat")
taliRain = as_tibble(taliRain$data) %>% dplyr::select(year = V1, dayIndex = V2, rain = V3)
taliRain %>% group_by(year) %>% summarise(a = sum(rain), sum(rain)) %>% summarise(mean(a))

taliRain = as_tibble(taliRain) %>% left_join(SynthRain  %>% filter(year<=500)%>% dplyr::select(myRain = rain, dayIndex, year, PET, - K, - SeriesNumber), by = c("year", "dayIndex"))

#aggregateReuslts---


test = ResultsTable %>% filter(isHolocene) %>% groupByParam %>% calculatePareto() # %>% filter(!pareto) %>% calculatePareto %>% filter(!pareto) %>% calculatePareto
 test %>% filter(pareto) %>% dplyr::select(sulfate, dustFlux) %>% gather() %>% ggplot(aes(y = value, x = key)) + geom_boxplot() + geom_point()+scale_y_continuous(breaks = scales::extended_breaks(20)) + theme(axis.title.x = element_blank())
test %>% filter(pareto) %>% dplyr::select(sulfate,dustFlux) %>% ungroup%>% summarise_all(median)

#Test real rain----
IMSRain$year = IMSRain$waterYear %% IMSRain$waterYear[1] + 1
IMSRainTest = IMSRain %>% left_join(SynthRain %>% dplyr::select(-rain) %>% filter(year <= max(IMSRain$year)), by = c("dayIndex", "year"))

pdf(file = paste("plots/", format(Sys.time(), "%b_%d_%Y_%H%M"), "ResultsSoilTallyTest.pdf"), width = 30, height = 16);
print(TAL);
dev.off()

gc(reset = TRUE)

#FC sensitivity ---
kiki1 = SynthRainE %>% filter(year == 1550, dayIndex >= 160 & dayIndex <= 161)
kiki1$rain[1] = 19
results1 = lapply(FCarray, FUN = function(X) CalcGypsum(kiki1, duration = 1, plotRes = 0, Depth = 50, AETFactor = 1, FieldCapacity = X, wieltingPoint = 0.014, thick = 2, verbose = 1));
kiki1$rain[1] = 4.2
results2 = lapply(FCarray, FUN = function(X) CalcGypsum(kiki1, duration = 1, plotRes = 0, Depth = 50, AETFactor = 0.6, FieldCapacity = X, wieltingPoint = 0.014, thick = 2, verbose = 1));
kiki1$rain[1] = 3.18
results5 = lapply(FCarray, FUN = function(X) CalcGypsum(kiki1, duration = 1, plotRes = 0, Depth = 50, AETFactor = 0.6, FieldCapacity = X, wieltingPoint = 0.014, thick = 2, verbose = 1));
kiki1$rain[1] = 3.79
results6 = lapply(FCarray, FUN = function(X) CalcGypsum(kiki1, duration = 1, plotRes = 0, Depth = 50, AETFactor = 0.6, FieldCapacity = X, wieltingPoint = 0.014, thick = 2, verbose = 1));
kiki1$rain[1] = 4.26
results7 = lapply(FCarray, FUN = function(X) CalcGypsum(kiki1, duration = 1, plotRes = 0, Depth = 50, AETFactor = 0.6, FieldCapacity = X, wieltingPoint = 0.014, thick = 2, verbose = 1));
kiki1$rain[1] = 3.55
results8 = lapply(FCarray, FUN = function(X) CalcGypsum(kiki1, duration = 1, plotRes = 0, Depth = 50, AETFactor = 0.6, FieldCapacity = X, wieltingPoint = 0.014, thick = 2, verbose = 1));


fctbl = tibble(FC = FCarray, EV1 = abs(results2 %>% map_dbl(~.x$maxWD) - 5.5), EV2 = abs(results5 %>% map_dbl(~.x$maxWD) - 4), EV3 = abs(results6 %>% map_dbl(~.x$maxWD) - 3.5), EV4 = abs(results7 %>% map_dbl(~.x$maxWD) - 4.5), EV5 = abs(results8 %>% map_dbl(~.x$maxWD) - 5)) %>% gather(group, rmsd, - FC) %>%
        mutate(rmsd = rmsd ^ 2) %>% group_by(FC) %>% summarise(rmsd = sqrt(sum(rmsd) / n()))
ggplot(fctbl, aes(x = FC, y = rmsd)) + geom_line() + geom_point() + coord_cartesian(ylim = c(0, 5), expand = 1) + ylab("RMSD [cm]") + xlab("FC [cm3/cm3]")
ggplot(tibble(rain = c(4.2, 3.18, 3.79, 4.26, 3.55), wd = c(5.5, 4, 3.5, 4.5, 5)), aes(rain, wd)) + geom_point()

#Kc calib----
kiki1 = SynthRainE %>% filter(year == 1550, dayIndex >= 160 & dayIndex <= 164)
kiki1$rain[1] = 4.6
kiki1$PET = 7.4
resultsKc = lapply(AETArray, FUN = function(X) CalcGypsum(kiki1, duration = 1, plotRes = 0, Depth = 50, AETFactor = X, FieldCapacity = 0.1, wieltingPoint = 0.013, thick = 10, verbose = 1));
fctbl = tibble(Kc = AETArray, ZEL11 = resultsKc %>% map_dbl(plotMoisture, c(2.5, 1.5, 1, 1, 1))) %>% gather(group, rmsd, - Kc)
ggplot(fctbl) + geom_line(aes(x = Kc, y = rmsd, color = group)) + scale_x_continuous(breaks = round(unique(fctbl$Kc), 3))


#sensitivity with gyp----
opt.AETF = 1.2;
opt.WP = 0.013;
opt.FC = 0.1;
opt.sulfate = 13;
opt.dust = 6;
seq.AETF = seq(0.8, 1.2, by = 0.1) * opt.AETF;
seq.WP = seq(0.8, 1.2, by = 0.1) * opt.WP;
seq.FC = seq(0.8, 1.2, by = 0.1) * opt.FC;
seq.rainSeq = seq(0.8, 1.2, by = 0.1) * opt.sulfate;
seq.dustSeq = seq(0.8, 1.2, by = 0.1) * opt.dust;
sens = as.matrix(bind_rows(crossing(AETF = seq.AETF, WP = opt.WP, FC = 0.1), crossing(AETF = 1.2, WP, FC = 0.1, rep = 1:15), crossing(AETF = 1.2, WP = 0.013, FC, rep = 1:15)))

results$T1.10 = c(results$T1.10, lapply(1:length(seq.dustSeq), function(X) CalcGypsum(SynthRainE %>% filter(year > 1000), duration = 13400, plotRes = 0, Depth = tail(T1.10Observed$bottom, 1), dustFlux = seq.dustSeq[X], random = F)))
results$T1.9 = c(results$T1.9, lapply(1:length(seq.dustSeq), function(X) CalcGypsum(SynthRainE %>% filter(year > 1000), duration = 11000, plotRes = 0, Depth = tail(T1.9Observed$bottom, 1), dustFlux = seq.dustSeq[X], random = F)))
results$zel11 = c(results$zel11, lapply(1:length(seq.dustSeq), function(X) CalcGypsum(SynthRainS %>% filter(year > 1000), duration = 10300, plotRes = 0, Depth = tail(Zel11Observed$bottom, 1), dustFlux = seq.dustSeq[X], random = F)))

results$T1.10 = c(results$T1.10, lapply(1:length(seq.FC), function(X) CalcGypsum(SynthRainE %>% filter(year > 1000), duration = 13400, plotRes = 0, Depth = tail(T1.10Observed$bottom, 1), FieldCapacity = seq.FC[X], random = F)))
results$T1.9 = c(results$T1.9, lapply(1:length(seq.FC), function(X) CalcGypsum(SynthRainE %>% filter(year > 1000), duration = 11000, plotRes = 0, Depth = tail(T1.9Observed$bottom, 1), FieldCapacity = seq.FC[X], random = F)))
results$zel11 = c(results$zel11, lapply(1:length(seq.FC), function(X) CalcGypsum(SynthRainS %>% filter(year > 1000), duration = 10300, plotRes = 0, Depth = tail(Zel11Observed$bottom, 1), FieldCapacity = seq.FC[X], random = F)))

results$T1.10 = c(results$T1.10, lapply(1:length(seq.AETF), function(X) CalcGypsum(SynthRainE %>% filter(year > 1000), duration = 13400, plotRes = 0, Depth = tail(T1.10Observed$bottom, 1), AETFactor = seq.AETF[X], random = F)))
results$T1.9 = c(results$T1.9, lapply(1:length(seq.AETF), function(X) CalcGypsum(SynthRainE %>% filter(year > 1000), duration = 11000, plotRes = 0, Depth = tail(T1.9Observed$bottom, 1), AETFactor = seq.AETF[X], random = F)))
results$zel11 = c(results$zel11, lapply(1:length(seq.AETF), function(X) CalcGypsum(SynthRainS %>% filter(year > 1000), duration = 10300, plotRes = 0, Depth = tail(Zel11Observed$bottom, 1), AETFactor = seq.AETF[X], random = F)))

results$T1.10 = c(results$T1.10, lapply(1:length(seq.WP), function(X) CalcGypsum(SynthRainE %>% filter(year > 1000), duration = 13400, plotRes = 0, Depth = tail(T1.10Observed$bottom, 1), wieltingPoint = seq.WP[X], random = F)))
results$T1.9 = c(results$T1.9, lapply(1:length(seq.WP), function(X) CalcGypsum(SynthRainE %>% filter(year > 1000), duration = 11000, plotRes = 0, Depth = tail(T1.9Observed$bottom, 1), wieltingPoint = seq.WP[X], random = F)))
results$zel11 = c(results$zel11, lapply(1:length(seq.WP), function(X) CalcGypsum(SynthRainS %>% filter(year > 1000), duration = 10300, plotRes = 0, Depth = tail(Zel11Observed$bottom, 1), wieltingPoint = seq.WP[X], random = F)))

results$T1.10 = c(results$T1.10, lapply(1:length(seq.rainSeq), function(X) CalcGypsum(SynthRainE %>% filter(year > 1000), duration = 13400, plotRes = 0, Depth = tail(T1.10Observed$bottom, 1), rainSO4 = seq.rainSeq[X], random = F)))
results$T1.9 = c(results$T1.9, lapply(1:length(seq.rainSeq), function(X) CalcGypsum(SynthRainE %>% filter(year > 1000), duration = 11000, plotRes = 0, Depth = tail(T1.9Observed$bottom, 1), rainSO4 = seq.rainSeq[X], random = F)))
results$zel11 = c(results$zel11, lapply(1:length(seq.rainSeq), function(X) CalcGypsum(SynthRainS %>% filter(year > 1000), duration = 10300, plotRes = 0, Depth = tail(Zel11Observed$bottom, 1), rainSO4 = seq.rainSeq[X], random = F)))

results$T1.10 = c(results$T1.10, lapply(1:length(seq.dustSeq), function(X) CalcGypsum(SynthRainE %>% filter(year > 1000), duration = 13400, plotRes = 0, Depth = tail(T1.10Observed$bottom, 1), dustFlux = seq.dustSeq[X], random = F)))
results$T1.9 = c(results$T1.9, lapply(1:length(seq.dustSeq), function(X) CalcGypsum(SynthRainE %>% filter(year > 1000), duration = 11000, plotRes = 0, Depth = tail(T1.9Observed$bottom, 1), dustFlux = seq.dustSeq[X], random = F)))
results$zel11 = c(results$zel11, lapply(1:length(seq.dustSeq), function(X) CalcGypsum(SynthRainS %>% filter(year > 1000), duration = 10300, plotRes = 0, Depth = tail(Zel11Observed$bottom, 1), dustFlux = seq.dustSeq[X], random = F)))
#create a table where all other parameter are fixed
#sens = ResultsTable %>% filter(isHolocene)
senstest = bind_rows(
WPSens = sens %>% filter(FC == opt.FC, AETF ==opt.AETF, sulfate == opt.sulfate, dustFlux == opt.dust, WP %in% seq.WP) %>% mutate(param = "\u03B8r", change = WP / opt.WP),
FCSens = sens %>% filter(WP == opt.WP, AETF == opt.AETF, sulfate == opt.sulfate, dustFlux == opt.dust, FC %in% seq.FC)  %>% mutate(param = "FC", change = FC / opt.FC),
AETFSens = sens %>% filter(FC == opt.FC, WP == opt.WP, sulfate == opt.sulfate, dustFlux == opt.dust, AETF %in% seq.AETF)  %>% mutate(param = "AET.F", change = AETF / opt.AETF),
sulfateSens = sens %>% filter(FC == opt.FC, AETF == opt.AETF, WP == opt.WP, dustFlux == opt.dust, sulfate %in% seq.rainSeq)   %>% mutate(param = "sulfate", change = sulfate / opt.sulfate),
dustFSens = sens %>% filter(FC == opt.FC, AETF == opt.AETF, sulfate == opt.sulfate, WP == opt.WP, dustFlux %in% seq.dustSeq)  %>% mutate(param = "dustFlux", change = dustFlux / opt.dust)
)



WP1 = senstest %>% filter(round(change, 1) %in% c(0.8, 1, 1.2)) %>% dplyr::select(change, Max_Concentration = PeakConc, Parameter = param) %>% gather("target", "value", - Parameter, - change) %>%
    ggplot(aes(y = value, x = factor(change), fill = Parameter)) + geom_bar(stat = "summary", fun = "mean", position = "dodge2") + facet_wrap(target ~ ., scales = "free_y", nrow = 3) +
    scale_y_continuous(name = "Concentration [meq/100g soil]") + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.text = element_text(size = 30), legend.title = element_text(size = 35), legend.position = c(), axis.title.y = element_text( hjust = 0.9)) + labs( color = "Parameter", x = "change from optimal") + scale_fill_brewer(type = "qual", palette = 7) + scale_x_discrete(labels = function(x) paste0(as.numeric(x) * 100, "%"))
# +geom_boxplot(outlier.shape = NA)
WP2 = senstest %>% filter(round(change, 1) %in% c(0.8, 1, 1.2)) %>% dplyr::select(change, Total_concentration = total, param) %>% gather("target", "value", - param, - change) %>%
    ggplot(aes(y = value, x = factor(change), fill = param)) + geom_bar(stat = "summary", fun = "mean", position = "dodge2") + facet_wrap(target ~ ., scales = "free_y", nrow = 3) +
    scale_y_continuous(name = "") + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.text = element_text(size = 30), legend.title = element_text(size = 35), legend.position = c("none")) + labs(color = "", x = "change from optimal") + scale_fill_brewer(type = "qual", palette = 7) + scale_x_discrete(labels = function(x) paste0(as.numeric(x) * 100, "%"))
# +geom_boxplot(outlier.shape = NA)
WP3 = senstest %>% filter(round(change, 1) %in% c(0.8, 1, 1.2)) %>% dplyr::select(change, Gypsum_depth = PeakDepth, param) %>% gather("target", "value", - param, - change) %>%
    ggplot(aes(y = value, x = factor(change), fill = param)) + geom_bar(stat = "summary", fun = "mean", position = "dodge2") + facet_wrap(target ~ ., scales = "free_y", nrow = 3) +
    scale_y_reverse(name = "Depth[cm]", position = "right") + theme(axis.text.x = element_text(size = 35), legend.text = element_text(size = 30), legend.title = element_text(size = 35), legend.position = c("none")) + labs(color = " change from \n optimal[%]", x = "Change from optimal") + scale_fill_brewer(type = "qual", palette = 7) + scale_x_discrete(labels = function(x) paste0(as.numeric(x) * 100, "%"))

    ggarrange(WP1,WP2,WP3)


rainDaysSens = resultsTable %>% filter(FC == opt.FC, AETF == opt.AETF, sulfate == opt.sulfate, WP == opt.WP, dustFlux == opt.dust, isHolocene, !is.na(rainDays))  %>%
     dplyr::select(rainDays, Gypsum_depth = PeakDepth, Total_concentration = total, Max_Concentration = PeakConc) %>% gather("target", "value", -rainDays) %>%
     ggplot(rainDaysSens, aes(x = rainDays, y = value)) +geom_point()+ facet_wrap(target ~ .)



RainArr = seq(12, 15, by = 1);
DustArr = c(seq(4,6, by = 0.5));
repetition = 1:10
rainDustArray = as.matrix(crossing(RainArr, DustArr, repetition))
nrow(rainDustArray)
#par computation----
source("evapostats/parallel.R");
RateTest = list()

resultsT1.10.par = lapply(1:nrow(rainDustArray), function(X) CalcGypsum(SynthRainE, duration = 13400, plotRes = 0, Depth = tail(T1.10Observed$bottom, 1), dustFlux = rainDustArray[X, 2], rainSO4 = rainDustArray[X, 1]))
resultsT1.9.par = lapply(1:nrow(rainDustArray), function(X) CalcGypsum(SynthRainE, duration = 11000, plotRes = 0, Depth = tail(T1.9Observed$bottom, 1), dustFlux = rainDustArray[X, 2], rainSO4 = rainDustArray[X, 1]));
resultsZel11.par = lapply(1:nrow(rainDustArray), function(X) CalcGypsum(SynthRainS, duration = 10300, plotRes = 0, Depth = tail(Zel11Observed$bottom, 1), dustFlux = rainDustArray[X, 2], rainSO4 = rainDustArray[X, 1]));


 lapply(1:nrow(rainDustArray), FUN = function(X) CalcGypsum(SynthRainE, duration = 100, plotRes = 0, Depth = 100, DustGyp=0.005,dustFlux = rainDustArray[X, 2], rainSO4 = rainDustArray[X, 1])))
save(RateTest, file =  "rateResults.Rdata")
bla = plotRate(RateTest)  + coord_cartesian(ylim = c())

 CVSHOnly= CV %>% filter(str_detect(profile, "T"));

dustTable = CV %>% filter(rain == 13 & profile == "T1.10" | rain == 14 & profile == "T1.9" | rain == 16 & profile == "zel11"  )
ggplot(dustTable, aes(x = dust, y = meanRMSD, color = profile, min = minRMSD, max = maxRMSD)) + geom_line() + geom_point() + geom_ribbon(alpha = 0.3) + geom_line(aes(y = normRMSD,size="- norm")) + coord_cartesian(ylim = c(), xlim = c())+
                labs(x = "dust flux [g/m2/year]", y = "normalised RMSD [-]", color = "profile", title = "dust flux using optimal rain sulfate ") 

rainTable = CV %>% filter(dust == 41& profile == "T1.10" | dust==33 & profile == "T1.9" | dust== 15 & profile == "zel11")
ggplot(rainTable, aes(x = rain, y = meanRMSD, group = profile, color = profile, max = maxRMSD, min = minRMSD)) + geom_line(size = 1) + geom_point() + geom_line(aes(y = normRMSD, size = "- norm")) + coord_cartesian(ylim = c(0,3),xlim = c()) + geom_ribbon(alpha = 0.3) +
    labs(x = "SO4 in rain [ml/L]", y = "normalised RMSD [-]", color = "profile", title = "sulfate in rain using optimal dust flux") # +stat_summary(aes(group=site,color = site),size = 2,geom = "line",fun.y = "mean") 

#cross validation
optimal = CV %>% filter(rain <20, dust < 10, n>1 ) %>% group_by(profile) %>% arrange(meanRMSD) %>% dplyr::slice(1)
optimal = LOOCV(CV %>% filter(dust < 10)) 
WP1 = ggplot(optimal, aes(x = rain, y = dust, color = profile)) + geom_point(size = 10) + labs(x = "sulfate [ml/l]", y = "dust flux [g/m2/yr]", fill = "RMSD [meq/100g soil]") + coord_cartesian(ylim = c(0, 50), xlim = c(10, 20))
cvres = LOOCV(CV %>% filter(dust < 10)) %>% transmute(All = val, profile) %>% gather("site", "value", - profile) # %>% left_join(LOOCV(CVSHOnly) %>% transmute(SH = val, profile), by = "profile") %>% left_join(LOOCV(CVZelOnly) %>% transmute(Zel = val, profile), by = "profile") 
WP2 = ggplot(cvres, aes(x = site, y = value)) + geom_point(stat = "summary", fun.y = "mean", size = 7) + geom_errorbar(stat = "summary", fun.data = "mean_se", fun.args = list(mult = 1.96))+geom_point(aes(color = profile),size = 5)+ labs(x = NULL, y = "norm. RMSD", fill = NULL,title = "CV results") #+ geom_bar(fill = "gray", stat = "identity", alpha = 0.1)
ggarrange(WP2, WP1,nrow=1)
#list of retrun variables (fucking genius)----
aWP1 = plotSoilResultsAgg(results$s10 %>% keep(~.x$RSO4 == 12.5 & .x$DF == 5.5), c(T1.10Observed %>% pull(mean)),  "T1.10 - Shehoret")
aWP2 = plotSoilResultsAgg(results$s9 %>% keep(~.x$RSO4 == 15 & .x$DF == 7), c(T1.9Observed %>% pull(mean)), "T1.9 - Shehoret")
aWP3 = plotSoilResultsAgg(results$zel11 %>% keep(~.x$RSO4 == optimal[3,1] & .x$DF == 6), c(Zel11Observed %>% pull(mean)), "zel11 - Zeelim")
aWP3 = plotSoilResultsAgg(bla, c(T1.1Observed[1:20,] %>% pull(mean)), "zel11 - Zeelim")
res = rbind(aWP1[[2]], aWP2[[2]], aWP3[[2]])
#WP4 = plotSoilResultsAgg(list(result), c(T2.1Observed  %>% pull(mean)), "T2.1")
#WP5 = plotSoilResultsAgg(results$zel13 %>% keep(~.x$RSO4 == 4 & .x$DF == 3.5), c(Zel13Observed %>% pull(mean)), "zel13")
plot_grid(aWP1, aWP2, aWP3, nrow = 2, legend = get_legend(aWP1),axis = "t")

WP1 = plotSoilResultsSurface(CV ,"zel11")
WP2 = plotSoilResultsSurface(CV , "T1.9")
WP3 = plotSoilResultsSurface(CV , "T1.10")
WP4 = plotSoilResultsSurface(CV ,c("zel11","T1.9","T1.10"))
res = rbind(WP1[[2]], WP2[[2]], WP3[[2]],WP4[[2]])


dustsens = CV %>% filter(round(dust,2) %in% c(5, 6, 7)&rain == 13) %>% group_by(par = dust) %>% summarise(value = joinRMSD(mean, n))
rainsens = CV %>% filter(round(rain,2) %in% c(10.4,13,15.6)&dust==6) %>% group_by(par=rain) %>% summarise(value = joinRMSD(mean, n))
diff(rainsens$value[1:2]) / diff(rainsens$par[1:2]) * rainsens$par[2]/rainsens$value[2] #* var(rainsens$par[1:2]) / var(rainsens$value[1:2])
diff(dustsens$value[1:2]) / diff(dustsens$par[1:2]) * dustsens$par[2] / dustsens$value[2] #* var(dustsens$par[1:2]) / var(dustsens$value[1:2])

diff(rainsens$value[3:2]) / diff(rainsens$par[3:2]) * rainsens$par[2] / rainsens$value[2] # * var(rainsens$par[2:3]) / var(rainsens$value[2:3])
diff(dustsens$value[2:3]) / diff(dustsens$par[2:3]) * dustsens$par[2] / dustsens$value[2] #* var(dustsens$par[2:3]) / var(dustsens$value[2:3])



for (item in unique(destsens$factor)) {
    temp = destsens %>% filter(factor == item)
    temp = (temp$value[2] - temp$value[1]) / (temp$rain[2] - temp$rain[1])
    resCV = resCV %>% add_row(name = item, val = temp)

}
#aetsENSITIVITY---
ggplot(tibble(meanWD = results$meanWD$value, Index30 = results$Index30$value, WDp80 = results$WDp80$value, FC = FCarray)) +
    geom_point(aes(FC, meanWD, color = "meanWD")) +
    geom_point(aes(FC, Index30, color = "Index30")) +
    geom_point(aes(FC, WDp80, color = "WDp80")) +
    scale_y_reverse() + ylab("depth")


ggplot(bla, aes(value, group = L1)) + stat_ecdf(aes(color = factor(L1))) + coord_flip(c(0,10)) + scale_x_reverse() + scale_y_reverse()



rm(list = ls(pattern = "(.par)"))
