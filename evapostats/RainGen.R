require(tidyverse)
require(tibble)
require(MASS)
require(actuar)
library(dplyr)
library(RODBC)
library(tictoc)

getDayAmount = function(Random, WetAfterDry, WetAfterWet, PrevDayAmount, RandomRainDepth) {
    getDayAmount = 0;
    WetProb = 0;
    #If prev was wet
    if (PrevDayAmount > 0) {
        Wetprob = WetAfterWet
    } else {
        WetProb = WetAfterDry
    }
    if (Random <= WetProb) {
        getDayAmount = RandomRainDepth
    }
    return(getDayAmount)
}

CalculateProbabilities = function() {
    #get wet after wet
    IMSRain$prevDay = IMSRain$dayIndex - lag(IMSRain$dayIndex, default = 0);
    IMSRain$WAW = IMSRain$prevDay == 1 | IMSRain$prevDay == -364;

    #histograms for wetdays
    WAWProb = hist(IMSRain %>% filter(WAW) %>% pull(dayIndex), right = FALSE, breaks = 1:366, plot = FALSE);
    WADProb = hist(IMSRain %>% filter(!WAW) %>% pull(dayIndex), right = FALSE, breaks = 1:366, plot = FALSE);
    WetProb = hist(IMSRain %>% pull(dayIndex), right = FALSE, breaks = 1:366, plot = FALSE);

    #Probability by dividing histograms
    PWAW = WAWProb$counts / (lag(WetProb$counts, default = tail(WetProb$counts, 1)));
    PWAD = WADProb$counts / (length(unique(IMSRain$year)) - lag(WetProb$counts, default = tail(WetProb$counts, 1)));

    #create probabilities for every annual day index
    Prob.Series = bind_cols("PWAW" = PWAW, "PWAD" = PWAD) %>% replace_na(list(PWAW = 0, PWAD = 0));

    return(Prob.Series)
}
GenerateSeries = function() {
    source("C:/Users/Lior/master/evapocalc/evapostats/PETGen.R")
    setwd('C:/Users/Lior/master/evapocalc/evapostats');

    con = odbcConnect("RainEvap")
    numOfyears = 1000;

    #get percipitation for 1988 to current
    IMSRain = tbl_df(sqlQuery(con, "SELECT *, year(time) as year,month(time) as month ,dayofyear(time) as dayIndex  FROM data_dream.precip_daily where (precip_daily.idstation = 347700 and year(precip_daily.time) >= 1958)"));

    #only rain days
    IMSRain = IMSRain[IMSRain$measure > 0.1,]

    pd = fitdistr(IMSRain$measure, "weibull");

    DaysProb = CalculateProbabilities();
       
    #rain generator ---

    DepthLimit = 0.1
    #create the synthetic rain series and pick random values for occurance and possible amount from weibull
    SynthRain = tibble(
                day = rep(1:365, numOfyears),
                year = rep(1:numOfyears, each = 365),
                RandomForOccurance = runif(length(day)),
                PotentialAmount = (rweibull(length(day), shape = pd$estimate[1], scale = pd$estimate[2])+DepthLimit),
                depth = 0);

    #First day, assume prev was dry
    SynthRain$depth[1] = getDayAmount(SynthRain$RandomForOccurance[1], DaysProb$PWAD[1], DaysProb$PWAW[1], 0, SynthRain$PotentialAmount[1])
    #iterating through the series based on previous row
    SynthRain$depth[2:nrow(SynthRain)] = sapply(2:nrow(SynthRain), function(X) getDayAmount(
                                           SynthRain$RandomForOccurance[X],
                                                DaysProb$PWAD[SynthRain$day[X]],
                                                DaysProb$PWAW[SynthRain$day[X]],
                                                SynthRain$depth[X-1],
                                                SynthRain$PotentialAmount[X]));

    finalRainSeries = PETGen(SynthRain);
    return(finalRainSeries)
    #---

}

plotResults = function() {
    #graphic---
    #bind syntetic rain and PET
    bla = SyntRain %>% group_by(year) %>% dplyr::summarise(annualRain = sum(depth),
                                                                annualPET = sum(PET));

    SynthRain$yeargroup = SynthRain$year %/% 30;

    chunks = SynthRain %>% group_by(yeargroup) %>% dplyr::summarise(std = sd(annualRain),
                                                                mean = mean(annualRain),
                                                                meanPET = mean(annualPET))


    blaIMS = IMSRain[which(IMSRain$year >= 1988),] %>% group_by(year) %>% dplyr::summarise(annualRain = sum(vals));
    IMSannualSD = sd(blaIMS$annualRain)
    IMSMeanAnnual = mean(blaIMS$annualRain)

    #bind syntetic rain and PET
    Synt = SyntRain %>% group_by(month);
    Synt = Synt %>% dplyr::summarise(
                                                     meanPET = mean(PET),
                                                     sumDepth = sum(depth),
                                                     stdPET = sd(PET)
                                                    )
    measured = RainSeries %>% group_by(month);
    measured = measured %>% dplyr::summarise(meanPET = mean(Pen),
                                                 stdPET = sd(Pen))


    #stdHistogram rain
    ggplot2::ggplot(data = bla, aes(std)) + geom_histogram(aes(y = ..density..)) + geom_density(aes(color = "Blue"), show.legend = FALSE) +
            geom_vline(aes(xintercept = IMSannualSD, color = "red")) + geom_vline(aes(xintercept = density(bla$std)$x[which.max(density(bla$std)$y)], color = "Blue"), linetype = "dashed") +
            labs(title = "STD histogram for 30 yr chunks\nMeasured STD is 13.6 ")
    mean(bla$std)

    #---
    IMSRain %>% group_by(year) %>% filter(year >1988 ) %>%summarise(sum = sum(measure)) %>% summarise(mean(sum))
    SynthRain %>% group_by(year) %>% summarise(annual = sum(depth)) %>% summarise(mean(annual))

    tali = as.tibble( read.csv("C:\\Users\\Lior\\master\\evapocalc\\DB\\RainSeriesEilatTest.csv"))
    tali %>% group_by(year) %>% summarise(annual = sum(depth)) %>% summarise(mean(annual))


    dailymean = IMSRain %>%  summarise(mean = mean(measure)) %>% pull(mean)
    meanEventsPerYear = IMSRain %>% group_by(year) %>% filter(year > 1988) %>% summarise(n = n()) %>% summarise(mean = mean(n)) %>% pull(mean)
    eventDurationMean = 1
}
