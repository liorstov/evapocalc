
load("C:/Users/liorst/Source/Repos/evapocalc/.RData")
require(tidyverse)
require(Rcpp)
require(ggplot2)
require(reshape2)
require(zoo)
require(tictoc)
require(R.matlab)
require(gganimate)


theme_set(theme_classic() + theme(legend.title = element_blank(), legend.key.size = unit(2, "picas"), legend.text = element_text(size = 15),
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


Rcpp::sourceCpp('C:/Users/liorst/source/repos/evapocalc/Calcyp/CSM.cpp', verbose = TRUE, rebuild = 0);
b <<- new(CSMCLASS);

duration = 1
result = CalcGypsum(SynthRainS %>% filter(year>5000), duration, plotRes = 0, Depth = 50, AETFactor =1, FieldCapacity = 0.1, wieltingPoint = 0.013, thick = 5, verbose = 1, dustFlux =6.52/10000 / 365, DustGyp = 0.005, rainCa = 35, rainSO4 =15);
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



FCarray = seq(0.02, 0.25, by = 0.01);
RainSO4Arr = seq(1, 20, length = 10);
wiltingPointArray = seq(0.01, 0.02,length = 10);
DustFluxArray = seq(from = 0.1,to =  3, length = 20);
AETArray = seq(from = 0.1, to = 2, length = 10);
RainFactorArray = seq(from = 0.05, to = 0.9, length = 2);
initIonArray = seq(from = 1, to = 20, length = 45);
AetRainComb =  as.matrix(crossing(RainFactorArray, AETArray))



#FC sensitivity ---
kiki1 = SynthRainE %>% filter(year == 1550, dayIndex >= 160 & dayIndex <= 161)
kiki1$rain[1] = 19
results1 = lapply(FCarray, FUN = function(X) CalcGypsum(kiki1, duration = 1, plotRes = 0, Depth = 50, AETFactor = 0.6, FieldCapacity = X, wieltingPoint = 0.014, thick = 2, verbose = 1, DustCa = 0, DustSO4 = 0, dustFlux = 0, rainCa = 0, rainSO4 = 0));
kiki1$rain[1] = 4.2
results2 = lapply(FCarray, FUN = function(X) CalcGypsum(kiki1, duration = 1, plotRes = 0, Depth = 50, AETFactor = 0.6, FieldCapacity = X, wieltingPoint = 0.014, thick = 2, verbose = 1, DustCa = 0, DustSO4 = 0, dustFlux = 0, rainCa = 0, rainSO4 = 0));
kiki1$rain[1] = 3.18
results5 = lapply(FCarray, FUN = function(X) CalcGypsum(kiki1, duration = 1, plotRes = 0, Depth = 50, AETFactor = 0.6, FieldCapacity = X, wieltingPoint = 0.014, thick = 2, verbose = 1, DustCa = 0, DustSO4 = 0, dustFlux = 0, rainCa = 0, rainSO4 = 0));
kiki1$rain[1] = 3.79
results6 = lapply(FCarray, FUN = function(X) CalcGypsum(kiki1, duration = 1, plotRes = 0, Depth = 50, AETFactor = 0.6, FieldCapacity = X, wieltingPoint = 0.014, thick = 2, verbose = 1, DustCa = 0, DustSO4 = 0, dustFlux = 0, rainCa = 0, rainSO4 = 0));
kiki1$rain[1] = 4.26
results7 = lapply(FCarray, FUN = function(X) CalcGypsum(kiki1, duration = 1, plotRes = 0, Depth = 50, AETFactor = 0.6, FieldCapacity = X, wieltingPoint = 0.014, thick = 2, verbose = 1, DustCa = 0, DustSO4 = 0, dustFlux = 0, rainCa = 0, rainSO4 = 0));
kiki1$rain[1] = 3.55
results8 = lapply(FCarray, FUN = function(X) CalcGypsum(kiki1, duration = 1, plotRes = 0, Depth = 50, AETFactor = 0.6, FieldCapacity = X, wieltingPoint = 0.014, thick = 2, verbose = 1, DustCa = 0, DustSO4 = 0, dustFlux = 0, rainCa = 0, rainSO4 = 0));


fctbl = tibble(FC = FCarray, EV1 = abs(results2 %>% map_dbl(~.x$maxWD) - 5.5), EV2 = abs(results5 %>% map_dbl(~.x$maxWD) - 4), EV3 = abs(results6 %>% map_dbl(~.x$maxWD) - 3.5), EV4 = abs(results7 %>% map_dbl(~.x$maxWD) - 4.5), EV5 = abs(results7 %>% map_dbl(~.x$maxWD) - 5)) %>% gather(group, rmsd, - FC) %>%
        mutate(rmsd = rmsd ^ 2) %>% group_by(FC) %>% summarise(rmsd = sqrt(sum(rmsd) / n()))
#ggplot(fctbl) + geom_line(aes(x = FC, y = rmsd, color = group)) + scale_x_continuous(breaks = round(unique(fctbl$FC), 3))+
ggplot(fctbl) + geom_line(aes(x = FC, y = rmsd))  


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
##for zel11
resultsSulfate = lapply(RainSO4Arr, FUN = function(X) CalcGypsum(SynthRainS, duration = 10000, plotRes = 0, Depth = 50, AETFactor = X, FieldCapacity = 0.1, wieltingPoint = 0.013, thick = 5, verbose = 0, dustFlux = 0.0006 / 365, rainCa = 35.58, rainSO4 = X));

#list of retrun variables (fucking genius)----
results = resultsSulfate %>% transpose %>% map_depth(2, ~ rowid_to_column(tibble(value = .x))) %>%
            pluck("gypsum") %>% melt(id.vars = "rowid", value.name = "value") %>% group_by(rowid) %>%
            dplyr::summarise(min = quantile(value, 0.05), gypsum = mean(value), max = quantile(value, 0.95)) %>% mutate(depth = (rowid - 0.5) * 5) %>%
            mutate(observed = c(3.141999936,10.38199997,14.75,14.75,14.75,14.75,14.75,14.75,14.75,14.74))
plotSoilResultsAgg(results)

#rmsd as compared to observed
rmsd = resultsSulfate %>% transpose %>% map_depth(2, ~ rowid_to_column(tibble(value = .x))) %>%
    pluck("gypsum") %>% map_dbl(~rmsd(.x, c(3.141999936, 10.38199997, 14.75, 14.75, 14.75, 14.75, 14.75, 14.75, 14.75, 14.74)))

ggplot(tibble(RainSO4Arr, rmsd), aes(x = RainSO4Arr, y = rmsd)) + geom_line() + geom_point() + scale_x_continuous(breaks = round(unique(RainSO4Arr), 3))


#aetsENSITIVITY---
ggplot(tibble(meanWD = results$meanWD$value, Index30 = results$Index30$value, WDp80 = results$WDp80$value, FC = FCarray)) +
    geom_point(aes(FC, meanWD, color = "meanWD")) +
    geom_point(aes(FC, Index30, color = "Index30")) +
    geom_point(aes(FC, WDp80, color = "WDp80")) +
    scale_y_reverse() + ylab("depth")


ggplot(bla, aes(value, group = L1)) + stat_ecdf(aes(color = factor(L1))) + coord_flip(c(0,10)) + scale_x_reverse() + scale_y_reverse()

#running FLUX AET combination
results = sapply(seq(1, nrow(AetRainComb)), FUN = function(X) CalcGypsum(years = 1000, RainFactor = AetRainComb[X, 1], AETFactor = AetRainComb[X, 2], observedArray = (Observed$zeelim.2EH)));
combi = results;
#running differrent inition
results = sapply(initIonArray, FUN = function(X) CalcGypsum(years = 10000, DustCa = 2, DustSO4 = 2, AETFactor = 90, InitialCa = X, initialSO4 = X, observedArray = (Observed$zeelim.2EH)));

#running differrent flux
results = sapply(DustFluxArray, FUN = function(X) CalcGypsum(years = 10000, DustCa = X, DustSO4 = X, observedArray = (Observed$zeelim.2EH)));


#running differrent WP
results = sapply(wiltingPointArray, FUN = function(X) CalcGypsum(years = 10000, DustCa = 2, DustSO4 = 2, AETFactor = 90, wieltingPoint = X, observedArray = (Observed$zeelim.2EH)));

#evapotranspiration
results = sapply(AETArray, FUN = function(X) CalcGypsum(years = 10000, AETFactor = X, Depth = 200, Getresults = FALSE, observedArray = (Observed$zeelim.2EH)));

#rain
results = sapply(RainFactorArray, FUN = function(X) CalcGypsum(years = 1000, RainFactor = X, Depth = 100, observedArray = (Observed$zeelim.2EH), raindata= EilatData$depth, PETData = EilatData$PET));

#tests
CalcGypsum(duration = 10000, RainFactor = 1, Depth = 200, observedArray = (Observed$shehoret1.MP), raindata = as.numeric(EilatData$depth), PETData = EilatData$PET, Getresults = TRUE, AETFactor = 1);

Observed$Calculated = obsCalc[, 2]
Observed[21:40,] = NA;

#convert to plotable
results = as.data.frame.array(t(results))
for (I in seq(1:ncol(results))) {
    results[, I] = unlist(results[I]);
}

ggplot(data = results, mapping = aes(x = results$DustCa, y = results$AETFactor, y = results$difference.observedArray.)) + geom_density();

#sesitivity test for InitIon
dustplot = ggplot(data = results) + geom_smooth(se = FALSE, mapping = aes(x = results$wieltingPoint, y = results$difference.observedArray.)) +
    labs(x = "Wilting Point [cm2]", y = "Gypsum accumulation depth difference [cm]",
    title = "sensitivity test for soil Initial soluble ions") + theme(text = element_text(size = 15));

#sesitivity test for WP
dustplot = ggplot(data = results) + geom_smooth(se = FALSE, mapping = aes(x = results$wieltingPoint, y = results$difference.observedArray.)) +
    labs(x = "Wilting Point [cm2]", y = "Gypsum accumulation depth difference [cm]",
    title = "sensitivity test for soil WiltingPoint") + theme(text = element_text(size = 15));


#sesitivity test for flux
dustplot = ggplot(data = results, mapping = aes(x = results$DustCa, y = results$difference.observedArray.)) + geom_point() + geom_path() +
   xlim (0,2.01)+
    labs(x = "yearly dust flux [g/m-2/yr-1]", y ="Gypsum accumulation depth difference [cm]",
    title = "sensitivity test for dust flux") + theme(text = element_text(size = 15, face = "bold"));

ggsave("plots/measuredOnly.png" )

#sensitivity test for AESFactor
ETPlot = ggplot(data = results, mapping = aes(x = results$AETFactor, y = results$difference.observedArray.)) + geom_point(se = FALSE) +
    labs(x = "Evapotranspiration Factor", y = "Gypsum accumulation depth difference [cm]",
    title = "sensitivity test for Evapotranspiration") + theme(text = element_text(size = 15, face = "bold")) + geom_path();

#sensitivity test for rain factor
ETPlot = ggplot(data = results, mapping = aes(x = results$RainFactor, y = results$difference.observedArray.)) + geom_point() +
    geom_point(mapping = aes(y = results$difference.observedArray.[7], x = results$RainFactor[7]), color = 'red', size = 2) +
    labs(x = "Rain Depth Factor", y = "Gypsum accumulation depth difference [cm]",
    title = "sensitivity test for Rain Factor \n AETFactor = 4.8; DustFlux = 1.5 g/m yr") +
    theme(text = element_text(size = 15), legend.text = element_text(size = 15, face = "bold")) + geom_path();


#response surfece
ggplot(results, aes(x = results$RainFactor, y = results$AETFactor, z = results$difference.observedArray)) +
    geom_raster(aes(fill = results$difference.observedArray), interpolate = TRUE, contour = TRUE) + geom_contour(color = "white") +    
    labs(x = "Rain Depth Factor", y = "Evaporation Factor", fill = "Gypsum accumulation\ndepth difference [cm]\n",
    title = "Response Surface for rain and evaporation factor") + theme(text = element_text(size = 15), legend.text = element_text(size = 12)) 

monthDF = as.data.frame(monthAET);
monthDF$month = factor(month.abb, levels = month.abb)
ggplot(data = monthDF, aes(x = month, y = monthAET)) + geom_bar(stat = "identity") +
    labs(y = "Average monthly AET [cm]", x = "" , title = "Average monthly AET for 10 ka" )

#plor gypsum profile
melted = melt(Observed, id.vars = "depth_roof")
melted$lineType = "solid";
melted$lineType[which(melted$variable == "Calculated")] = "dashed";
ggplot(data = melted[which(melted$variable %in% c( "shehoret1.MP", "shehoret3.MP", "zeelim.12H", "zeelim.13H", "zeelim.11MH" ,"zeelim.2EH", "zeelim.1EH")),],
    mapping = aes(x = value, y = depth_roof,, group = variable, colour = variable, linetype = lineType)) + geom_path() + scale_y_reverse() +
    labs(x = "Gypsum [meq/100 gr soil]", y = "Depth [cm]", colour = "Site name", title = "Measured Data") +
    scale_linetype_manual(values = c("solid", "solid")) + guides(linetype = FALSE) + theme(text = element_text(size = 20));


#"zeelim.12H", "zeelim.13H", "zeelim.11MH" ,"zeelim.2EH", "zeelim.1EH", "Calculated"

#find minimum value
minimal = temp[3, temp[1,] == min(unlist(temp[1,]))]

##create for Eilta
EilatData = read.csv("Eilat.csv",);

##get amount values 
EilatAmountCdf = ecdf(SedomData[,3]);

##Interval 
#EilatInterval = tail(EilatData$time, -1) - head(EilatData$time, -1);
#EilatInterval = EilatInterval[which(EilatInterval < 180)];

##get interval CDF
#EilatIntervalCDF = ecdf(EilatInterval);
#Intervalcdf = cbind(get("x", environment(EilatIntervalCDF)), get("y", environment(EilatIntervalCDF)));

##get CDF function
Amountcdf = cbind(get("x", environment(EilatAmountCdf)), get("y", environment(EilatAmountCdf)));
test = as.data.frame.array(Amountcdf);
ggplot(test, aes(x = test$V1)) + stat_ecdf(geom = "point") + labs(x = "Event magnitude [mm]", title = "ecdf for Sodom station") + theme(text = element_text(size = 30))

test


##generate Rain series

#RainSeries = numeric(NumberOfYears * 365);
#interval = 0;
#RainSeries = sapply(RainSeries, FUN = function(X) GetRandomValue(X));

##create Module from CSM
#Rcpp::sourceCpp('C:/Users/Lior/master/evapocalc/Calcyp/main.cpp');
#n = new(CSM,2);


