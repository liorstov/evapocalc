
load("C:/Users/liorst/Source/Repos/evapocalc/.RData")
require(tidyverse)
require(Rcpp)
require(ggplot2)
require(reshape2)
require(zoo)
require(tictoc)
require(R.matlab)
require(gganimate)
require(pbapply)


theme_set(theme_classic() + theme( legend.key.size = unit(2, "picas"), legend.text = element_text(size = 15),
axis.text.x = element_text(size = 20, angle = 43,hjust = 1),
axis.text.y = element_text(size = 25),
axis.title.y = element_text(size = 20),
axis.title.x = element_text(size = 20),
panel.grid.major = element_line(colour = "grey93"),
panel.grid.minor = element_line(colour = "grey93"),
panel.border = element_rect(colour = "black", fill = NA, size = 0)))
Sys.setlocale("LC_TIME", "English_Israel.1255");

#set wd 
setwd("C:/Users/liorst/source/repos/evapocalc");

source("evapostats/Functions.R");
source("evapostats/PETGen.R");
source("evapostats/RainGen.R");
Rcpp::sourceCpp('C:/Users/liorst/source/repos/evapocalc/Calcyp/CSM.cpp', verbose = TRUE, rebuild = 0);
b <<- new(CSMCLASS);

#stationElat = 347700;
#stationElatEvap = 347704;
#stationSedom = 337000;
IMSRain = GetImsRain(station = 347700, stationEvap = 347704);
rainSeriesResults = GenerateSeries(NumOfSeries = 1000, IMSRain = IMSRain);
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



duration = 1000
result = CalcGypsum(SynthRainS, duration, plotRes = 1, Depth = 50, AETFactor =1, FieldCapacity = 0.1, wieltingPoint = 0.013, thick = 5, verbose = 0, dustFlux =6.52/10000 / 365, DustGyp = 0.005, rainCa = 35, rainSO4 =15);
plotSoilResults(result)


#sandbox

kiki1 = SynthRainE %>% filter(year == 1550, dayIndex >= 160 & dayIndex <= 161)


taliRain = readMat("DB\\RainTali.mat")
taliRain = as_tibble(taliRain$data) %>% dplyr::select(year = V1, dayIndex = V2, rain = V3)
taliRain %>% group_by(year) %>% summarise(a = sum(rain), sum(rain)) %>% summarise(mean(a))

taliRain = as_tibble(taliRain) %>% left_join(SynthRain  %>% filter(year<=500)%>% dplyr::select(myRain = rain, dayIndex, year, PET, - K, - SeriesNumber), by = c("year", "dayIndex"))

#Test real rain
IMSRain$year = IMSRain$waterYear %% IMSRain$waterYear[1] + 1
IMSRainTest = IMSRain %>% left_join(SynthRain %>% dplyr::select(-rain) %>% filter(year <= max(IMSRain$year)), by = c("dayIndex", "year"))

pdf(file = paste("plots/", format(Sys.time(), "%b_%d_%Y_%H%M"), "ResultsSoilTallyTest.pdf"), width = 30, height = 16);
print(TAL);
dev.off()

gc(reset = TRUE)

#----






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
#ggplot(fctbl) + geom_line(aes(x = FC, y = rmsd, color = group)) + scale_x_continuous(breaks = round(unique(fctbl$FC), 3))+
ggplot(fctbl, aes(x = FC, y = rmsd)) + geom_line() + geom_point() + coord_cartesian(ylim = c(0, 5), expand = 1) + ylab("RMSD [cm]") + xlab("FC [cm3/cm3]")


#Kc calib----
kiki1 = SynthRainE %>% filter(year == 1550, dayIndex >= 160 & dayIndex <= 164)
kiki1$rain[1] = 4.6
kiki1$PET =7.4
resultsKc = lapply(AETArray, FUN = function(X) CalcGypsum(kiki1, duration = 1, plotRes = 0, Depth = 50, AETFactor = X, FieldCapacity = 0.1, wieltingPoint = 0.013, thick = 10, verbose = 1));
fctbl = tibble(Kc = AETArray, ZEL11 = resultsKc %>% map_dbl(plotMoisture, c(2.5, 1.5, 1, 1, 1))) %>% gather(group, rmsd, - Kc)
ggplot(fctbl) + geom_line(aes(x = Kc, y = rmsd, color = group)) + scale_x_continuous(breaks = round(unique(fctbl$Kc), 3))



test = tibble(SMeanWD = results %>% map(~pluck(., c("SMeanWD"))) %>% unlist(),
                Index3 = results %>% map(~pluck(., c("Index03"))) %>% unlist(),
                SWDp80 = results %>% map(~pluck(., c("SWDp80"))) %>% unlist()) %>% rowid_to_column() %>% gather(Index, depth, - rowid)
ggplot(test) + geom_bar(aes(x  = rowid * 70, y = depth, color = Index)) + scale_y_reverse(expand = c(0, 0.0), limits = c(50, 0)) +
                        scale_x_continuous(expand = c(0, 0.0), name = "years")


#Rain Sulfate---
observedProfiles = read_csv("DB\\Calcyp_GypsoilHorizons_caso4.slab.csv") %>% dplyr::select(1:12) %>% drop_na()
Zel11Observed = observedProfiles %>% filter(str_detect(SiteName, "11")) 
Zel12Observed = observedProfiles %>% filter(str_detect(SiteName, "12")) 
Zel13Observed = observedProfiles %>% filter(str_detect(SiteName, "13")) 
T1.9Observed = observedProfiles %>% filter(str_detect(SiteName, "T1-9"))
T1.10Observed = observedProfiles %>% filter(str_detect(SiteName, "T1-10"))

FCarray = seq(0.02, 0.3, by = 0.01);
RainArr = c(seq(0, 11.96, length = 5), seq(11.97, 14.63, length = 10), seq(14.64, 25, length = 5));
DustArr = c(seq(0, 1.16, length = 5), seq(1.17,1.43, length = 10), seq(1.44, 7, length = 5));
rainDustArray = as.matrix(crossing(RainArr, DustArr))

##for zel11

resultsDustT1.10 = pblapply(1:nrow(rainDustArray), FUN = function(X) CalcGypsum(SynthRainE, duration = 13400, plotRes = 0, Depth = tail(T1.10Observed$bottom, 1), dustFlux = rainDustArray[X, 2], rainSO4 = rainDustArray[X, 1]));
resultsDustT1.9 = pblapply(1:nrow(rainDustArray), FUN = function(X) CalcGypsum(SynthRainE, duration = 13400, plotRes = 0, Depth = tail(T1.9Observed$bottom, 1), dustFlux = rainDustArray[X, 2], rainSO4 = rainDustArray[X, 1]));
resultsRainZel11 = pblapply(1:nrow(rainDustArray), FUN = function(X) CalcGypsum(SynthRainS, duration = 10300, plotRes = 0, Depth = tail(Zel11Observed$bottom, 1), dustFlux = rainDustArray[X, 2], rainSO4 = rainDustArray[X, 1]));
resultsDustZel12 = pblapply(1:nrow(rainDustArray), FUN = function(X) CalcGypsum(SynthRainS, duration = 8900, plotRes = 0, Depth = tail(Zel12Observed$bottom, 1), dustFlux = rainDustArray[X, 2], rainSO4 = rainDustArray[X, 1]));
resultsDustZel13 = pblapply(1:nrow(rainDustArray), FUN = function(X) CalcGypsum(SynthRainS, duration = 2800, plotRes = 0, Depth = tail(Zel13Observed$bottom, 1), dustFlux = rainDustArray[X, 2], rainSO4 = rainDustArray[X, 1]));
save.image(file = "test.RData")
#list of retrun variables (fucking genius)----
selectRes =  as_tibble(rainDustArray) %>% rowid_to_column %>% filter(RainArr == RMSDRain.T.10$rain[67]) %>% pull(rowid)
plotSoilResultsAgg(resultsDustT1.10[67], c(T1.10Observed %>% pull(mean)))
plotSoilResultsAgg(resultsDustT1.9[74], c(T1.9Observed %>% pull(mean)))
plotSoilResultsAgg(resultsRainZel11[77], c(Zel11Observed %>% pull(mean)))
plotSoilResultsAgg(resultsDustZel12[31], c(Zel12Observed %>% pull(mean)))

RMSDRain.T.10 = plotSoilResultsMean(resultsDustT1.10, c(T1.10Observed %>% pull(mean)), rainDustArray) %>% rename(T1.10 = value)
RMSDRain.T1.9 = plotSoilResultsMean(resultsDustT1.9, c(T1.9Observed %>% pull(mean)), rainDustArray) %>% rename(T1.9 = value)
RMSDRain.Zel11 = plotSoilResultsMean(resultsRainZel11, c(Zel11Observed %>% pull(mean)), rainDustArray) %>% rename(zel11 = value)
RMSDRain.Zel12 = plotSoilResultsMean(resultsDustZel12, c(Zel12Observed %>% pull(mean)), rainDustArray) %>% rename(zel12 = value)
RMSDRain.Zel13 = plotSoilResultsMean(resultsDustZel13, c(Zel13Observed %>% pull(mean)), rainDustArray) %>% rename(zel13 = value)
#rmsd as compared to observed
RMSDRain.T.10 = plotSoilResultsMean(resultsRainT.10, c(T1.10Observed %>% pull(mean)), RainArr) 
RMSDRain.T1.9 = plotSoilResultsMean(resultsRainT1.9, c(T1.9Observed %>% pull(mean)), RainArr)
RMSDRain.Zel11 = plotSoilResultsMean(resultsRainZel11, c(Zel11Observed %>% pull(mean)), RainArr)
RMSDRain.Zel12 = plotSoilResultsMean(resultsRainZel12, c(Zel12Observed %>% pull(mean)), RainArr)
bla = plotSoilResultsRMSD(resultsRainZel11Const, c(Zel11Observed %>% pull(mean)), RainArr)
bla = plotSoilResultsRMSD(resultsDustZel12, c(Zel12Observed %>% pull(mean)), DustArr)
bla = plotSoilResultsRMSD(resultsRainZel12Const, c(Zel12Observed %>% pull(mean)), RainArr)
bla = plotSoilResultsRMSD(resultsDustZel11, c(Zel11Observed %>% pull(mean)), DustArr)
bla = plotSoilResultsRMSD(resultsDustZel11Const, c(Zel11Observed %>% pull(mean)), DustArr)
bla = plotSoilResultsRMSD(resultsDustZel12Const, c(Zel12Observed %>% pull(mean)), DustArr)
bla = plotSoilResultsRMSD(resultsSulfateT1.10, c(T1.10Observed %>% pull(mean)), DustArr)
RMSDRain.T.10[which.min(bla1$),]

#combine results of rmsd
bla1 = RMSDRain.T.10 %>%    filter(dust == RMSDRain.T.10$dust[67]) %>% dplyr::select(-dust, T1.10 = value) #0.74
bla2 = RMSDRain.T1.9 %>%    filter(dust == RMSDRain.T1.9$dust[74]) %>% dplyr::select(-dust, T1.9 = value) #5.08
bla3 = RMSDRain.Zel11 %>%   filter(dust == RMSDRain.Zel11$dust[77]) %>% dplyr::select(-dust, Zel11 = value) #2.59
bla4 = RMSDRain.Zel12 %>%   filter(dust == RMSDRain.Zel12$dust[31]) %>% dplyr::select(-dust, Zel12 = value) #5.04
bla5 = RMSDRain.Zel13 %>%   filter(dust == RMSDRain.Zel13$dust[43]) %>% dplyr::select(-dust, Zel13 = value) #5.04

rmsdTable = bla1 %>% left_join(bla2, by = "rain") %>% left_join(bla3, by = "rain") %>% left_join(bla4, by = "rain") %>% left_join(bla5, by = "rain") %>% mutate(zel = sqrt((Zel11 ^ 2 + Zel12 ^ 2 + Zel13 ^ 2) / 3), SH = sqrt((T1.9 ^ 2 + T1.10 ^ 2) / 2)) %>% gather("factor", "value", - rain) # 
ggplot(rmsdTable %>% filter(factor %in% c("zel","SH")), aes(x = rain, y = value, color = factor)) + geom_line() + geom_point() +
                labs(x = "SO4 in rain [ml/L]", y = "RMSD meq/100g soil", color = "profile", title = "sulfate in rain using optimal dust flux")
ggplot(rmsdTable, aes(x = rain, y = value, color = factor)) + geom_line() + geom_point() +
                labs(x = "SO4 in rain [ml/L]", y = "RMSD meq/100g soil", color = "profile", title = "sulfate in rain using optimal dust flux")

#cross validation
CV = RMSDRain.T.10 %>% left_join(RMSDRain.T1.9, by = c("rain", "dust")) %>% left_join(RMSDRain.Zel11, by = c("rain", "dust")) %>% left_join(RMSDRain.Zel12, by = c("rain", "dust")) %>% left_join(RMSDRain.Zel13, by = c("rain", "dust")) %>%
            gather("factor", "value", - rain, - dust) %>% unnest() %>% mutate(RMSD = sqrt((value) / (n)))

resCV = tibble(name = numeric(), val = numeric(), dust = numeric(), rain = numeric())
for (item in unique(CV$factor)) {
    temp = CV %>% filter(factor != item) %>% group_by(rain, dust) %>% summarise(RMSD = sqrt(sum(value) / sum(n)))
    temp = temp %>% arrange(RMSD) %>% head(1) 

    res = CV %>% filter(factor == item, dust == temp$dust, rain == temp$rain)
    resCV = resCV %>% add_row(name = item, val = res$RMSD, dust = res$dust,rain = res$rain)

}
    {

        

}
CV %>% group_by(rain, dust) %>% summarise(sqrt(sum(value) / sum(n)))


#i = 0;
#for (item in 3:ncol(CV)) { i = i+1
    #headline = CV %>% mutate(s = sqrt(rowSums(.[, c(-1, -2, - as.integer(item))]) / (CVn %>% mutate(den = rowSums(.[, c(-1, -2, - as.integer(item))])) %>% pull(den)))) %>%
            #arrange(s) %>% head(1)
    #CVRes[i] = (sqrt(headline %>% pull(item) / CVn[1, item]))
#}

blak2 = RMSDRain.T1.9 %>% filter(rain == RMSDRain.T1.9$rain[74]) %>% dplyr::select(-rain, T1.9 = value) #5.08
blak3 = RMSDRain.Zel11 %>% filter(rain == RMSDRain.Zel11$rain[77]) %>% dplyr::select(-rain, Zel11 = value) #2.59
blak4 = RMSDRain.Zel12 %>% filter(rain == RMSDRain.Zel12$rain[31]) %>% dplyr::select(-rain, Zel12 = value) #5.04

rmsdTable = blak1 %>% left_join(blak2, by = "dust") %>% left_join(blak3, by = "dust") %>% left_join(blak4, by = "dust") %>% gather("factor", "value", - dust) #%>% mutate(value = sqrt((value.x ^ 2 + value.y ^ 2) / 2)) 
ggplot(rmsdTable, aes(x = dust, y = value, color = factor)) + geom_line() + geom_point() +
                labs(x = "dust flux [g/m2/year]", y = "RMSD meq/100g soil", color = "profile", title = "dust flux usin optimal rain")


ggplot(blak1, aes(x = dust, y = T1.10)) + geom_line() + geom_point() +    labs(x = "dust flux [g/m2/year]", y = "RMSD meq/100g soil")



#aetsENSITIVITY---
ggplot(tibble(meanWD = results$meanWD$value, Index30 = results$Index30$value, WDp80 = results$WDp80$value, FC = FCarray)) +
    geom_point(aes(FC, meanWD, color = "meanWD")) +
    geom_point(aes(FC, Index30, color = "Index30")) +
    geom_point(aes(FC, WDp80, color = "WDp80")) +
    scale_y_reverse() + ylab("depth")


ggplot(bla, aes(value, group = L1)) + stat_ecdf(aes(color = factor(L1))) + coord_flip(c(0,10)) + scale_x_reverse() + scale_y_reverse()

