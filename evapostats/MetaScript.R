
#load("C:/Users/liorst/Source/Repos/evapocalc/.RData")
load("C:/Users/liorst/Source/Repos/evapocalc/synth.RData")

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


opt.AETF <<- 1.2;
opt.WP <<- 0.013;
opt.FC <<- 0.1;
opt.sulfate <<- 10;
opt.dust <<- 2;
opt.WHC <<- 0.087
seq.WHC = seq(0.8, 1.2, by = 0.2) * opt.WHC %>% rep(60);
seq.AETF = seq(0.8, 1.2, by = 0.2) * opt.AETF %>% rep(60);
seq.WP = seq(0.8, 1.2, by = 0.2) * opt.WP %>% rep(60);
seq.FC = seq(0.8, 1.2, by = 0.2) * opt.FC %>% rep(60);
seq.rainSeq = seq(0.8, 1.2, by = 0.2) * opt.sulfate %>% rep(60);
seq.dustSeq = seq(0.8, 1.2, by = 0.2) * opt.dust %>% rep(60);

theme_set(theme_classic() + theme(legend.key.size = unit(1, "line"), legend.text = element_text(size = 25),
axis.text.x = element_text(size = 28, angle = 43, hjust = 1), title = element_text(size = 20),
strip.text = element_text( size = 20),
axis.text.y = element_text(size = 28),
axis.title.y = element_text(size = 35,margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.text.y.right = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
axis.title.x = element_text(size = 35),
legend.key.height = unit(3, "line"),
panel.grid.major = element_line(colour = "grey93"),
panel.grid.minor = element_line(colour = "grey93"),
panel.border = element_rect(colour = "black", fill = NA, size = 0))) + guides(color = guide_legend(override.aes = list(size = 8)))
Sys.setlocale("LC_TIME", "English_Israel.1255");

#set wd 
setwd("C:/Users/liorst/source/repos/evapocalc");


source("evapostats/Functions.R");
source("evapostats/PETGen.R");
source("evapostats/RainGen.R");
loadObsProfiles()

Rcpp::sourceCpp('C:/Users/liorst/source/repos/evapocalc/Calcyp/CSM.cpp', verbose = TRUE, rebuild = 0);
cppModule <- new(CSMCLASS);

#rain gen ----
#stationElat = 347700;
#stationElatEvap = 347704;
#stationSedom = 337000;
factor = tibble(key = numeric(), value = numeric())
for (item in seq(0.1,1,by = 0.1)) {
    test = GenerateSeries(station = 347700, stationEvap = 347704, NumOfSeries = 60, PETfactor = item)
    value = test %>% group_by(year) %>% summarise_all(sum) %>% summarise_all(mean) %>% pull(PET)
    factor = factor %>% add_row(key = item, value = value)
}
factorS = tibble(key = numeric(), value = numeric())
for (item in seq(0.1, 1, by = 0.1)) {
    test = GenerateSeries(station = 337000, stationEvap = 337000, PETfactor = 1, NumOfSeries = 60)
    value = test %>% group_by(year) %>% summarise_all(sum) %>% summarise_all(mean) %>% pull(PET)
    factorS = factorS %>% add_row(key = item, value = value)
}


SynthRainEP= GenerateSeries(station = 347700, stationEvap = 347704, NumOfSeries = 1000,  AnuualRain = 70, WetDays = 15)
IMSRain = GetImsRain(station = 337000, stationEvap = 337000)
results = plotResults(SynthRain, IMSRain, rainSeriesResults$DaysProb, PETresults$PETProb, 1);

IMSRain %>% filter(rain > 0) %>% group_by(waterYear) %>% summarise(s = sum(rain), n = n(),PET = sum(pen)) %>% drop_na() %>% summarise(mean(n), mean(s), mean(PET,na.rm = T))

result = CalcGypsum(SynthRainE, SynthRainEP, duration = 60000, Depth = 150, rainstat = F, random = T, withRunoff = T, withFC = T,rainSO4 = 13)
plotSoilResults(results[[1]], T1.1,T)


#tests---
test = lapply(1:100, FUN = function(X) CalcGypsum(SynthRainE, SynthRainEP, duration = 62000-5000, Depth = 150, rainstat = F, random = T))

testOSL1 = parLapply(cl, 1:1, fun = function(X) CalcGypsum(SynthRainE, SynthRainEP, duration = 62500, Depth = 150, rainstat = F, random = T, withFC = T, withRunoff = T, rainSO4 = 55)) 
testOSL3 = parLapply(cl, 1:40, fun = function(X) CalcGypsum(SynthRainE, SynthRainEP, duration = (62000 + 5000), Depth = 150, rainstat = F, random = T, withFC = F, withRunoff = T))
testOSL2 = parLapply(cl, 1:40, fun = function(X) CalcGypsum(SynthRainE, SynthRainEP, duration = (62000 - 5000), Depth = 150, rainstat = F, random = T, withFC = F, withRunoff = T))

test1 = CalcGypsum(SynthRainS, SynthRainS, duration = 60000, Depth = 150, rainstat = T, random = F, withRunoff = T, withFC = T, getWD = T)
test2 = CalcGypsum(SynthRainS, SynthRainS, duration = 60000, Depth = 150, rainstat = T, random = F, withRunoff = T, withFC = F, getWD = T)
test3 = CalcGypsum(SynthRainS, SynthRainS, duration = 60000, Depth = 150, rainstat = T, random = F, withRunoff = F, withFC = T, getWD = T)
test4 = CalcGypsum(SynthRainS, SynthRainS, duration = 60000, Depth = 150, rainstat = T, random = F, withRunoff = F, withFC = F, getWD = T)
test1$module = "FC and runoff"
 test2$module = "runoff only"
 test3$module = "FC only"
 test4$module = "no FC no runoff"
plotSoilResults(test1, Zel1)

rect = c(testOSL1, testOSL2, testOSL3) %>% RectanglingResults(T1.1)
rect %>% group_by(duration) %>% summarise(min(total), max(total), min(PeakDepth), max(PeakDepth))
rect = rect %>% rowwise() %>% mutate(xmin = 40, xmax = 60.5, ymin = 47.5, ymax = 47.5, inrange = (between(total, xmin, xmax) & between(PeakDepth, ymin, ymax))) %>% ungroup()
rect %>% group_by(duration) %>% summarise(sum(as.integer(inrange))/ n(), mean(total), median(PeakDepth))
rect %>% ggplot(aes(factor(duration), total)) + geom_boxplot(outlier.color = NA) + coord_cartesian() + scale_y_reverse() #+ geom_rect(aes(xmin = 40, xmax = 60.5, ymin = 47.5, ymax = 47.5), fill = NA, color = "black")

#print Data of measured profiles
observedProfiles %>% group_by(SiteName) %>% summarise(gyp = sumGypsum(caso4, 5), depth = GetGypsicHorizonDepth(tibble(caso4, top, bottom)), age = mean(AvgAge)) %>% ggplot(aes(age, depth, color = SiteName, shape = grepl("Sh", SiteName))) + geom_point(size = 8)

kiki1 = SynthRainE %>% filter(year == 1550, dayIndex >= 160 & dayIndex <= 161)
    

taliRain = readMat("DB\\RainTali.mat")
taliRain = as_tibble(taliRain$data) %>% dplyr::select(year = V1, dayIndex = V2, rain = V3)
taliRain %>% group_by(year) %>% summarise(a = sum(rain), sum(rain)) %>% summarise(mean(a))

taliRain = as_tibble(taliRain) %>% left_join(SynthRain  %>% filter(year<=500)%>% dplyr::select(myRain = rain, dayIndex, year, PET, - K, - SeriesNumber), by = c("year", "dayIndex"))

#aggregateReuslts---

test = resultsTable %>% filter(isHolocene) %>% groupByParam %>% calculatePareto() # %>% filter(!pareto) %>% calculatePareto %>% filter(!pareto) %>% calculatePareto
 bla %>% filter(pareto) %>% dplyr::select(sulfate, dustFlux) %>% gather() %>% ggplot(aes(y = value, x = key)) + geom_boxplot() + geom_point()+scale_y_continuous(breaks = scales::extended_breaks(20)) + theme(axis.title.x = element_blank())
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
seq.AETF = seq(0.8, 1.2, by = 0.1) * opt.AETF %>% rep(150);
seq.WP = seq(0.8, 1.2, by = 0.1) * opt.WP %>% rep(150);
seq.FC = seq(0.8, 1.2, by = 0.1) * opt.FC %>% rep(150);
seq.rainSeq = seq(0.8, 1.2, by = 0.1) * opt.sulfate %>% rep(150);
seq.dustSeq = seq(0.8, 1.2, by = 0.1) * opt.dust %>% rep(150);
sens = as.matrix(bind_rows(crossing(AETF = seq.AETF, WP = opt.WP, FC = 0.1), crossing(AETF = 1.2, WP, FC = 0.1, rep = 1:15), crossing(AETF = 1.2, WP = 0.013, FC, rep = 1:15)))

for (i in 1:30) {


results$T1.10 = c(results$T1.10, lapply(1:length(seq.dustSeq), function(X) CalcGypsum(SynthRainE %>% filter(year > 1000), duration = 13400, plotRes = 0, Depth = tail(T1.10Observed$bottom, 1), dustFlux = seq.dustSeq[X])))
results$T1.9 = c(results$T1.9, lapply(1:length(seq.dustSeq), function(X) CalcGypsum(SynthRainE %>% filter(year > 1000), duration = 11000, plotRes = 0, Depth = tail(T1.9Observed$bottom, 1), dustFlux = seq.dustSeq[X])))
results$zel11 = c(results$zel11, lapply(1:length(seq.dustSeq), function(X) CalcGypsum(SynthRainS %>% filter(year > 1000), duration = 10300, plotRes = 0, Depth = tail(Zel11Observed$bottom, 1), dustFlux = seq.dustSeq[X])))

results$T1.10 = c(results$T1.10, lapply(1:length(seq.FC), function(X) CalcGypsum(SynthRainE %>% filter(year > 1000), duration = 13400, plotRes = 0, Depth = tail(T1.10Observed$bottom, 1), FieldCapacity = seq.FC[X])))
results$T1.9 = c(results$T1.9, lapply(1:length(seq.FC), function(X) CalcGypsum(SynthRainE %>% filter(year > 1000), duration = 11000, plotRes = 0, Depth = tail(T1.9Observed$bottom, 1), FieldCapacity = seq.FC[X])))
results$zel11 = c(results$zel11, lapply(1:length(seq.FC), function(X) CalcGypsum(SynthRainS %>% filter(year > 1000), duration = 10300, plotRes = 0, Depth = tail(Zel11Observed$bottom, 1), FieldCapacity = seq.FC[X])))

results$T1.10 = c(results$T1.10, lapply(1:length(seq.AETF), function(X) CalcGypsum(SynthRainE %>% filter(year > 1000), duration = 13400, plotRes = 0, Depth = tail(T1.10Observed$bottom, 1), AETFactor = seq.AETF[X])))
results$T1.9 = c(results$T1.9, lapply(1:length(seq.AETF), function(X) CalcGypsum(SynthRainE %>% filter(year > 1000), duration = 11000, plotRes = 0, Depth = tail(T1.9Observed$bottom, 1), AETFactor = seq.AETF[X])))
results$zel11 = c(results$zel11, lapply(1:length(seq.AETF), function(X) CalcGypsum(SynthRainS %>% filter(year > 1000), duration = 10300, plotRes = 0, Depth = tail(Zel11Observed$bottom, 1), AETFactor = seq.AETF[X])))

results$T1.10 = c(results$T1.10, lapply(1:length(seq.WP), function(X) CalcGypsum(SynthRainE %>% filter(year > 1000), duration = 13400, plotRes = 0, Depth = tail(T1.10Observed$bottom, 1), wieltingPoint = seq.WP[X])))
results$T1.9 = c(results$T1.9, lapply(1:length(seq.WP), function(X) CalcGypsum(SynthRainE %>% filter(year > 1000), duration = 11000, plotRes = 0, Depth = tail(T1.9Observed$bottom, 1), wieltingPoint = seq.WP[X])))
results$zel11 = c(results$zel11, lapply(1:length(seq.WP), function(X) CalcGypsum(SynthRainS %>% filter(year > 1000), duration = 10300, plotRes = 0, Depth = tail(Zel11Observed$bottom, 1), wieltingPoint = seq.WP[X])))

results$T1.10 = c(results$T1.10, lapply(1:length(seq.rainSeq), function(X) CalcGypsum(SynthRainE %>% filter(year > 1000), duration = 13400, plotRes = 0, Depth = tail(T1.10Observed$bottom, 1), rainSO4 = seq.rainSeq[X])))
results$T1.9 = c(results$T1.9, lapply(1:length(seq.rainSeq), function(X) CalcGypsum(SynthRainE %>% filter(year > 1000), duration = 11000, plotRes = 0, Depth = tail(T1.9Observed$bottom, 1), rainSO4 = seq.rainSeq[X])))
results$zel11 = c(results$zel11, lapply(1:length(seq.rainSeq), function(X) CalcGypsum(SynthRainS %>% filter(year > 1000), duration = 10300, plotRes = 0, Depth = tail(Zel11Observed$bottom, 1), rainSO4 = seq.rainSeq[X])))

}
#create a table where all other parameter are fixed
#sens = ResultsTable %>% filter(isHolocene)

senstest = bind_rows(
WPSens = resultsTable %>% filter(FC == opt.FC, AETF ==opt.AETF, sulfate == opt.sulfate, dustFlux == opt.dust, WP %in% seq.WP) %>% mutate(param = "\u03B8r", change = WP / opt.WP),
FCSens = resultsTable %>% filter(WP == opt.WP, AETF == opt.AETF, sulfate == opt.sulfate, FC %in% seq.FC) %>% mutate(param = "FC", change = FC / opt.FC),
AETFSens = resultsTable %>% filter(FC == opt.FC, WP == opt.WP, sulfate == opt.sulfate, dustFlux == opt.dust, AETF %in% seq.AETF) %>% mutate(param = "AET.F", change = AETF / opt.AETF),
sulfateSens = resultsTable %>% filter(FC == opt.FC, AETF == opt.AETF, WP == opt.WP, dustFlux == opt.dust, sulfate %in% seq.rainSeq) %>% mutate(param = "sulfate", change = sulfate / opt.sulfate),
dustFSens = resultsTable %>% filter(FC == opt.FC, AETF == opt.AETF, sulfate == opt.sulfate, WP == opt.WP, dustFlux %in% seq.dustSeq) %>% mutate(param = "dustFlux", change = dustFlux / opt.dust)
)
bla = senstest %>% rowid_to_column(var = "id")%>% filter(change == 1) %>% group_by(PeakConc) %>% summarise(mean(FC), n = n(), target = mean(id)) %>% filter(n == 1) %>% pull(target)
senstest = senstest %>% slice(-bla)

senstest = senstest %>% filter(change == 1) %>% dplyr::select(refMaxconc = PeakConc, refDepth = PeakDepth, refTotal = total) %>% summarise_all(median) %>% cbind(senstest) %>% tibble
senstest = senstest %>% mutate(refMaxconc = PeakConc / refMaxconc, refDepth = PeakDepth / refDepth, refTotal = total / refTotal)
save(senstest,file = "senstable.RData")
WP1 = senstest %>% filter(round(change, 1) %in% c(0.8, 1, 1.2)) %>% dplyr::select(change, Max_Concentration = PeakConc, Parameter = param) %>% gather("target", "value", - Parameter, - change) %>%
    ggplot(aes(y = value, x = factor(change), fill = Parameter)) + stat_boxplot(geom = 'errorbar', linetype = 1, width = 0.5, position = position_dodge(0.75)) + geom_boxplot() + geom_point(stat = "summary", fun = "mean", shape = 5, size = 4, position = position_dodge(0.75)) + facet_wrap(target ~ ., scales = "free_y", nrow = 3) +
    scale_y_continuous(name = "Concentration [meq/100g soil]") + theme( strip.text = element_text(size = 30), legend.text = element_text(size = 30), legend.title = element_text(size = 35), axis.title.y = element_text( hjust = 0.9)) + labs( color = "Parameter", x = "change from optimal") + scale_fill_brewer(type = "qual", palette = 7) + scale_x_discrete(labels = function(x) paste0(as.numeric(x) * 100, "%"))
# +geom_boxplot(outlier.shape = NA)
WP2 = senstest %>% filter(round(change, 1) %in% c(0.8, 1, 1.2)) %>% dplyr::select(change, Total_concentration = total, param) %>% gather("target", "value", - param, - change) %>%
    ggplot(aes(y = value, x = factor(change), fill = param)) + stat_boxplot(geom = 'errorbar', linetype = 1, width = 0.5, position = position_dodge(0.75)) + geom_boxplot() + geom_point(stat = "summary", fun = "mean", shape = 5, size = 4, position = position_dodge(0.75)) + facet_wrap(target ~ ., scales = "free_y", nrow = 3) +
    scale_y_continuous(name = "") + theme(strip.text = element_text(size = 30), legend.text = element_text(size = 30), legend.title = element_text(size = 35)) + labs(color = "", x = "change from optimal") + scale_fill_brewer(type = "qual", palette = 7) + scale_x_discrete(labels = function(x) paste0(as.numeric(x) * 100, "%"))
# +geom_boxplot(outlier.shape = NA)
WP3 = senstest %>% filter(round(change, 1) %in% c(0.8, 1, 1.2)) %>% dplyr::select(change, Gypsum_depth = PeakDepth, param) %>% gather("target", "value", - param, - change) %>%
    ggplot(aes(y = value, x = factor(change), fill = param)) + stat_boxplot(geom = 'errorbar', linetype = 1, width = 0.5, position = position_dodge(0.75)) + geom_boxplot() + geom_point(stat = "summary", fun = "mean", shape = 5, size = 4, position = position_dodge(0.75)) + facet_wrap(target ~ ., scales = "free_y", nrow = 3) +
    scale_y_reverse(name = "Depth[cm]", position = "right") + theme(strip.text = element_text(size = 30), axis.text.x = element_text(size = 35), legend.text = element_text(size = 30), legend.title = element_text(size = 35)) + labs(color = " change from \n optimal[%]", x = "Change from optimal") + scale_fill_brewer(type = "qual", palette = 7) + scale_x_discrete(labels = function(x) paste0(as.numeric(x) * 100, "%"))

ggarrange(WP1 + theme(axis.text.x = element_blank(), axis.title.x = element_blank(),legend.key.size =  unit(5,"line")),
    WP2 + theme(axis.text.x = element_blank(), axis.title.x = element_blank(),legend.position = c("")),
        WP3 + theme( legend.position = c("")))



rainDaysTable = resultsTable %>% filter(profile == "T1.9", FC == opt.FC, AETF == opt.AETF, sulfate == opt.sulfate, dustFlux == opt.dust, WP == opt.WP, !is.na(rainDays)) %>%
    group_by(meanDay = round(AnnualRain / rainDays,1)) %>% summarise(median = median(total), min = quantile(total, 0.2), max = quantile(total, 0.8)) %>% ggplot(aes(x = meanDay, y = median, ymax = max, ymin = min)) + geom_line() + geom_ribbon(alpha = 0.5)



## sens test for weather series---
sensTest = resultsTable %>% filter(dustFlux == 11, FC == opt.FC, AETF == opt.AETF, WP == opt.WP, sulfate == 12) %>% group_by(AnnualRain = round(AnnualRain), rainDays = round(rainDays), region) %>% summarise_if(is.numeric,median,) %>% mutate(ObsRain = ifelse(region == "shehoret", 27.5, 41), ObsDepth = ifelse(region == "shehoret", 27.5, 22.5), ObsConc = ifelse(region == "shehoret", 2.5, 10), ObsTotal = ifelse(region == "shehoret", 2.5, 10), meanDay = (AnnualRain / rainDays), Gypsum_depth = PeakDepth, Max_Concentration = PeakConc, Total_concentration = total) %>%
    gather("target", "value", Gypsum_depth , Max_Concentration, Total_concentration )
sensTest %>% filter(target == "Gypsum_depth") %>% ggplot(aes(x = AnnualRain, y = value, color = q99, size = rainDays)) + geom_point(shape = 1) + facet_grid(region ~ target) + scale_color_gradientn(colours = rainbow(5, rev = T)) + labs(x = "Annual rain [mm]", y = "depth [cm]", color = "q99 [mm]", size = "#  annual\nrain days") + scale_y_reverse() + geom_point(aes(x = ObsRain, y = ObsDepth), color = "black", size = 10)

sensTest %>% filter(target != "Gypsum_depth")  %>% ggplot(aes(x = AnnualRain, y = value, color = (rainDays))) + geom_point() + facet_wrap(region ~ target, scales = "free_y") + scale_color_gradientn(colours = rainbow(5, rev = T)) + labs(x = "Annual rain [mm]", y = "gypsum[meq/100g soil]", color = "q99 [mm]", size = "#  annual\nrain days")
#---


RainArr = seq(20, 25, by = 1);
DustArr = c(seq(5,25, by = 5));
repetition = 1:5
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
#WP4 = plotSoilResultsAgg(list(result), c(T2.1  %>% pull(mean)), "T2.1")
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
