
load("C:/Users/liorst/Source/Repos/evapocalc/.RData")
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


theme_set(theme_classic() + theme(legend.key.size = unit(2, "picas"), legend.text = element_text(size = 25), legend.title = element_text(size = 25),
axis.text.x = element_text(size = 28, angle = 43, hjust = 1), title = element_text(size = 20),
axis.text.y = element_text(size = 28),
axis.title.y = element_text(size = 35),
axis.title.x = element_text(size = 35),
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
plotSoilResults(resultsDustT1.10[[1]])



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
RainArr = seq(10, 15, by = 1);
DustArr = c(seq(6,100, by = 1));
repetition = 1:6
rainDustArray = as.matrix(crossing(RainArr, DustArr, repetition))
nrow(rainDustArray)
#par computation----
source("evapostats/parallel.R");


 CVSHOnly= CV %>% filter(str_detect(profile, "T"));

dustTable = CV %>% filter(rain == 13 & profile == "T1.10" | rain == 14 & profile == "T1.9" | rain == 16 & profile == "zel11"  )
ggplot(dustTable, aes(x = dust, y = meanRMSD, color = profile, min = minRMSD, max = maxRMSD)) + geom_line() + geom_point() + geom_ribbon(alpha = 0.3) + geom_line(aes(y = normRMSD,size="- norm")) + coord_cartesian(ylim = c(), xlim = c())+
                labs(x = "dust flux [g/m2/year]", y = "normalised RMSD [-]", color = "profile", title = "dust flux using optimal rain sulfate ") 

rainTable = CV %>% filter(dust == 41& profile == "T1.10" | dust==33 & profile == "T1.9" | dust== 15 & profile == "zel11")
ggplot(rainTable, aes(x = rain, y = meanRMSD, group = profile, color = profile, max = maxRMSD, min = minRMSD)) + geom_line(size = 1) + geom_point() + geom_line(aes(y = normRMSD, size = "- norm")) + coord_cartesian(ylim = c(0,3),xlim = c()) + geom_ribbon(alpha = 0.3) +
    labs(x = "SO4 in rain [ml/L]", y = "normalised RMSD [-]", color = "profile", title = "sulfate in rain using optimal dust flux") # +stat_summary(aes(group=site,color = site),size = 2,geom = "line",fun.y = "mean") 

#cross validation
optimal = CV %>% filter(dust < 10) %>% group_by(profile) %>% arrange(meanRMSD) %>% slice(1)
optimal = LOOCV(CV %>% filter(dust < 96)) 
WP1 = ggplot(optimal, aes(x = rain, y = dust, color = profile)) + geom_point(size = 6, show.legend = FALSE) + labs(x = "sulfate [ml/l]", y = "dust flux [g/m2/yr]", title = "validation parameters", fill = "RMSD [meq/100g soil]") + coord_cartesian(ylim = c(0, 90), xlim = c(0, 20))
cvres = LOOCV(CV %>% filter(dust < 10)) %>% transmute(All = val, profile) %>% gather("site", "value", - profile) # %>% left_join(LOOCV(CVSHOnly) %>% transmute(SH = val, profile), by = "profile") %>% left_join(LOOCV(CVZelOnly) %>% transmute(Zel = val, profile), by = "profile") 
WP2 = ggplot(cvres, aes(x = site, y = value)) + geom_point(stat = "summary", fun.y = "mean", size = 7) + geom_errorbar(stat = "summary", fun.data = "mean_se", fun.args = list(mult = 1.96))+geom_point(aes(color = profile),size = 5)+ labs(x = NULL, y = "norm. RMSD", fill = NULL,title = "CV results") #+ geom_bar(fill = "gray", stat = "identity", alpha = 0.1)
ggarrange(WP2, WP1,nrow=1)
#list of retrun variables (fucking genius)----
aWP1 = plotSoilResultsAgg(results$s10 %>% keep(~.x$RSO4 == 11 & .x$DF == 41), c(T1.10Observed %>% pull(mean)),  "T1.10")
aWP2 = plotSoilResultsAgg(results$s9 %>% keep(~.x$RSO4 == 14 & .x$DF == 33), c(T1.9Observed %>% pull(mean)), "T1.9")
aWP3 = plotSoilResultsAgg(results$zel11 %>% keep(~.x$RSO4 == 16& .x$DF == 15), c(Zel11Observed %>% pull(mean)), "zel11")
#WP4 = plotSoilResultsAgg(results$zel12 %>% keep(~.x$RSO4 == 5 & .x$DF == 5.5), c(Zel12Observed %>% pull(mean)), "zel12")
#WP5 = plotSoilResultsAgg(results$zel13 %>% keep(~.x$RSO4 == 4 & .x$DF == 3.5), c(Zel13Observed %>% pull(mean)), "zel13")
ggarrange(aWP1,aWP2,aWP3,WP4,nrow = 2)



#i = 0;
#for (item in 3:ncol(CV)) { i = i+1
    #headline = CV %>% mutate(s = sqrt(rowSums(.[, c(-1, -2, - as.integer(item))]) / (CVn %>% mutate(den = rowSums(.[, c(-1, -2, - as.integer(item))])) %>% pull(den)))) %>%
            #arrange(s) %>% head(1)
    #CVRes[i] = (sqrt(headline %>% pull(item) / CVn[1, item]))
#}

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
