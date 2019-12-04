
require(tidyverse)
require(Rcpp)
require(ggplot2)
require(reshape2)
require(zoo)


theme_set(theme_classic() + theme(legend.title = element_blank(), legend.key.size = unit(2, "picas"), legend.text = element_text(size = 15),
axis.text.x = element_text(size = 20, angle = 43,hjust = 1),
axis.text.y = element_text(size = 25),
axis.title.y = element_text(size = 20),
axis.title.x = element_text(size = 20),
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
PETresults = PETGen(rainSeriesResults$SynthRain, IMSRain,40);

SynthRain = rainSeriesResults$SynthRain;
SynthRain$PET = PETresults$SynthPET;
SynthRain$K = PETresults$K;
PETProb = PETresults$PETProb;
rainProb = rainSeriesResults$DaysProb;
results = plotResults(SynthRain, IMSRain, rainSeriesResults$DaysProb, PETresults$PETProb, 1);

pdf(file = paste("plots/", format(Sys.time(), "%b_%d_%Y_%H%M"), "Results.pdf"), width = 30, height = 16);
print(results);
dev.off()

Observed = as.data.frame(read.csv("DB/measured.CSV"));
Observed = RawData2Compartments(Observed, 5);
Observed$depth_roof = Observed$compartment * 5;


temp = SedomData
newval = (temp[which(temp[, 3] > 20, 3), 3]) * 0.3
temp[which(temp[, 1] > 20, 1), 1] = newval

Sys.setenv(PATH = "%PATH%;C:/Rtools/gcc-4.6.3/bin;c:/Rtools/bin")
Rcpp::sourceCpp('C:/Users/liorst/source/repos/evapocalc/Calcyp/CSM.cpp', verbose = TRUE);
wiltingPointArray = seq(0.001, 0.1,length = 45);
DustFluxArray = seq(from = 0.1,to =  2, length = 20);
AETArray = seq(from = 1, to = 5, length = 20);
RainFactorArray = seq(from = 0.05, to = 0.9, length = 2);
initIonArray = seq(from = 1, to = 20, length = 45);
AetRainComb =  as.matrix(crossing(RainFactorArray, AETArray))

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
bla = CalcGypsum(duration = 20, raindata = SynthRain$rain, PETData = SynthRain$PET * 0, DustCa = 0.5, DustSO4 = 0.5, RainFactor = 1, Depth = 100, Getresults = TRUE);

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


