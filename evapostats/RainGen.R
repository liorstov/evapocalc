require(tidyverse)
require(tibble)
require(MASS)
require(actuar)
library(dplyr)
library(RODBC)
library(tictoc)
require(fitdistrplus)
require(egg)
library(lubridate)
require(pracma)


getDayAmount = function(dayIndex, Random, WetAfterDry, WetAfterWet, PrevDayAmount) {
    getDayAmount = 0;
    WetProb = 0;
    #If prev was wet
    if (dayIndex == 1) {
        return(0);
    }
    if (PrevDayAmount > 0.1) {
        Wetprob = WetAfterWet
    } else {
        WetProb = WetAfterDry
    }
    if (Random < WetProb) {
        getDayAmount = 0;# rweibull(1, shape = 0.4851, scale = 1.3458) + 0.1;
    }
    return(getDayAmount)
}

#returns tibble with station parameters
getStationParam = function(isStat) {
    weibParams = tibble(stationName = "Eilat", stationNumber = 347700, scale = 1.3458, shape = 0.4851) %>%
    add_row(stationName = "Sedom", stationNumber= 337000, shape = 0.5994, scale = 1.61)

    return(weibParams %>% filter(stationNumber == isStat))
}

CalculateProbabilities = function(IMSRain) {
    #get wet after wet
    IMSRainOnly = IMSRain %>% filter(rain > 0.1)
    IMSRainOnly$prevDay = IMSRainOnly$dayIndex - lag(IMSRainOnly$dayIndex, default = 0);
    IMSRainOnly$WAW = IMSRainOnly$prevDay == 1 | IMSRainOnly$prevDay == -364;

    #histograms for wetdays
    WAWProb = hist(IMSRainOnly %>% filter(WAW) %>% pull(dayIndex), right = FALSE, breaks = 1:366, plot = FALSE);
    WADProb = hist(IMSRainOnly %>% filter(!WAW) %>% pull(dayIndex), right = FALSE, breaks = 1:366, plot = FALSE);
    WetProb = hist(IMSRainOnly %>% pull(dayIndex), right = FALSE, breaks = 1:366, plot = FALSE);

    #Probability by dividing histograms
    PWAW = WAWProb$counts / (lag(WetProb$counts, default = tail(WetProb$counts, 1)));
    PWAD = WADProb$counts / (length(unique(IMSRainOnly$waterYear)) - lag(WetProb$counts, default = tail(WetProb$counts, 1)));
    PWET = WetProb$counts / length(unique(IMSRainOnly$waterYear));

    #create probabilities for every annual day index
    Prob.Series = bind_cols("PWETOrig" = PWET, "PWAWOrig" = PWAW, "PWADOrig" = PWAD, dayIndex = 1:365) %>% replace_na(list(PWAWOrig = 0, PWADOrig = 0));
    Prob.Series = Prob.Series %>% add_column(month = month(dmy("1-09-2000") + Prob.Series$dayIndex));
    #smoothing
    window1 <- 50
    window2 = 40  
    Prob.Series = Prob.Series %>%
        mutate(PWET = stats::filter(stats::filter(Prob.Series$PWETOrig, rep(1 / window1, window1), method = "con", sides = 2, circular = 1), rep(1 / window2, window2), method = "con", circular = 1, sides = 2),
               PWAW = stats::filter(stats::filter(Prob.Series$PWAWOrig, rep(1 / window1, window1), method = "con", sides = 2, circular = 1), rep(1 / window2, window2), method = "con", circular = 1, sides = 2),
               PWAD = stats::filter(stats::filter(Prob.Series$PWADOrig, rep(1 / window1, window1), method = "con", sides = 2, circular = 1), rep(1 / window2, window2), method = "con", circular = 1, sides = 2))

   
    return(Prob.Series)
}
GetImsRain = function(station = 337000) {
    #stationElat = 347700;
    #stationSedom = 337000;
    con = odbcConnect("RainEvap")
    IMSRain = tbl_df(sqlQuery(con, paste("SELECT *,  year(time) as year,month(time) as month, measure as rain FROM data_dream.precip_daily where ((precip_daily.idstation = 337000))")));
    EvapSeries = tbl_df(sqlQuery(con, "SELECT *, measure as pen, year(time) as year,month(time) as month ,dayofyear(time) as dayIndex FROM data_dream.pet_daily where idstation = 337000 and isnull(measure) = 0"));
    RODBC::odbcCloseAll();

    #expend to get full year representation
    IMSRain = IMSRain %>% full_join(tibble(time = seq(IMSRain$time[1], IMSRain$time[1]+365, by = '1 day')), by = "time")

    if (nrow(EvapSeries)) {
        IMSRain = IMSRain %>% full_join(EvapSeries, by = "time") %>% dplyr::select(time, rain, pen);
    }
    IMSRain$rain = IMSRain %>% pull(rain) %>% replace_na(replace = 0);

    #creating day Index of water year and removing leap years
    IMSRain = IMSRain %>% mutate(waterYear = ifelse(month(time) %in% 9:12, year(time) + 1, year(time))) %>%
                                 mutate(dayIndex = as.numeric(difftime(time, make_date(waterYear-1,9), units = "days"))+1) %>% filter(dayIndex <= 365)
    measuredYears <<- length(unique(IMSRain$waterYear))
    #IMSRain = IMSRain %>% filter(idstation == station)

    #IMSRain = IMSRain %>% filter(year >= 1948 ,year <= 2017, measure > 0.1)
    

    return(IMSRain)
}
GenerateSeries = function(NumOfSeries = 1000, withEvapo = FALSE, IMSRain)
{
    DaysProb = CalculateProbabilities(IMSRain);    
    DepthLimit = 0.1
    NumOfSeries <<- NumOfSeries;
    #weibull parameters are taken for francesco
    station = getStationParam(IMSRain$idstation[1])

    
    numOfyears <<- measuredYears * NumOfSeries
    #create the synthetic rain series and pick random values for occurance and possible amount from weibull
    #rain generator ---
    tic()

    #create matrix with weibull values for depth
    weibull = matrix(rweibull(365 * numOfyears, station$shape, station$scale) + DepthLimit, nrow = numOfyears, ncol = 365);

    #matrix with random values 
    randMat = matrix(runif(365 * numOfyears), nrow = numOfyears, ncol = 365);
    SynthRain = matrix(0, nrow = numOfyears, ncol = 365, dimnames = list(1:numOfyears, 1:365))
    
    for (days in 2:ncol(SynthRain)) {
        SynthRain[, days] = as.numeric((SynthRain[, days - 1] & (randMat[, days] < DaysProb$PWAW[days])) | (!SynthRain[, days - 1] & (randMat[, days] < DaysProb$PWAD[days])));
    }

    #addd amounts
    SynthRain[which(SynthRain == 1)] = weibull[which(SynthRain == 1)];
    toc()

    #as matrix
    RainMat <<- SynthRain;

    SynthRain = as_tibble(melt(SynthRain, value.name = "rain", varnames = c("year", "dayIndex")))
    #clear(c("weibull", "randMat"));

    if (withEvapo) {
        tic()
        print("calc PET")
        SynthRain$PET = PETGen(SynthRain, IMSRain);
        print(paste("PET: ", toc()))
    }

    return(SynthRain) 
    #---
}

Hist2tibble = function(histogram) {
    histogram = histogram %>% modify_at(1, ~ tail(.x, -1))
    return (tibble(breaks = histogram$breaks, counts = histogram$counts, density = histogram$density))
}

plotResults = function(SynthRain, IMSRain, withEvapo = FALSE) {
    #graphic---   
    #dividing the sim series to grups 
    SynthRain$SeriesNumber = SynthRain$SeriesNumber = rep(1:NumOfSeries, each = measuredYears)
    SynthRain$Wet = (SynthRain$rain > 0.1) * 1;
    #wet day prob----
   

    SimRainDay = SynthRain %>% dplyr::select(year, dayIndex, Wet) %>% spread(dayIndex, Wet, drop = FALSE) 

    SimRainDay = as_tibble(melt(SimRainDay[, 2:ncol(SimRainDay)])) %>% group_by(day = variable) %>% 
                                dplyr::summarise(min = quantile(value, 0.05), median = mean(value), max = quantile(value, 0.95)) %>%
                                    add_column(avg = (SynthRain %>% group_by(dayIndex) %>% dplyr::summarise(prob = (sum(Wet) / numOfyears)) %>% pull(prob))) %>%
                                   add_column(measure = IMSRain %>% group_by(day = dayIndex) %>% mutate(rain = 1 * (rain > 0.1)) %>% dplyr::summarise(Measured = (sum(rain) / measuredYears)) %>% pull(Measured)) %>%
                                    add_column(month = (dmy("1-09-2000") + 0:364))

    p1 = ggplot(data = SimRainDay, aes(x = month, group = 1)) + geom_line(aes(y = measure, color = "Measured")) + ylab("Wet Probability [-]") +
                                geom_line(aes(y = avg, color = "Calculated"), size = 1.5) + ggtitle(paste(numOfyears, "years")) +
                               # geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.2) +
                                scale_x_date(date_labels = "%B")
    #----
     
    #annual rain----

    #histogram for observed
    tic()
    plotLimit = 120;
    histBreaksSize = 5;
    IMSRainAnn = IMSRain %>% filter(rain > 0.1) %>% group_by(waterYear) %>% dplyr::summarise(annual = sum(rain), WetDays = n()) %>% filter(annual <= plotLimit)
    IMShistannual = hist(IMSRainAnn$annual, breaks = seq(0, plotLimit, histBreaksSize), plot = 0) %>%
                        Hist2tibble() %>% mutate(obsDensity = counts / measuredYears)
    #histogram for observed
    simRainAnn = SynthRain %>% filter(rain > 0.1) %>% group_by(year, SeriesNumber) %>% dplyr::summarise(annual = sum(rain), WetDays = n()) %>% filter(annual <= plotLimit)
    rainHist = simRainAnn %>% group_by(SeriesNumber) %>%
                       dplyr::summarise(a = list(hist(annual, breaks = seq(0, plotLimit, histBreaksSize), plot = 0)), n = n()) %>% pull(a) %>%
                        map(~Hist2tibble(.x)) %>% reduce(bind_rows)

    #aggregating 
    k = rainHist %>% group_by(breaks) %>% dplyr::summarise(simulated = mean(density), min = quantile(density, 0.05), median = quantile(density, 0.5), max = quantile(density, 0.95)) %>%
                left_join(IMShistannual, by = "breaks")

    p2 = ggplot(k, aes(x = breaks)) + geom_line(aes(y = simulated, color = "Calculated"), size = 1.5) + xlim(0, plotLimit) +
                                geom_line(data = IMShistannual, aes(y = density, color = "Measured")) +
                                geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.5) + xlab("Annual rain mm/year")
    toc()

   
    #---
    #annual wet days
    p3 = ggplot2::ggplot(data = simRainAnn, aes(WetDays)) + stat_density(aes(color = "Calculated"), size = 1.5, geom = "line", bw = 0.5) +
         stat_density(data = IMSRainAnn, aes(WetDays, color = "Measured"), geom = "line", bw = 0.5) + scale_x_continuous(breaks = seq(0, 30, 5)) + xlab("# Annual wet days")
    #---
    #PET per month
    p4=p5=ggplot()
    if (withEvapo) {


        DailyPET = IMSRain %>% filter(!is.na(pen)) %>% group_by(dayIndex) %>% dplyr::summarise(Measured = mean(pen))
        PETquant = SynthRain %>% group_by(dayIndex) %>% dplyr::summarise(median = mean(PET, 0.5))
        SimPet = SynthRain %>% group_by(dayIndex) %>% dplyr::summarise(PET = mean(PET)) %>% full_join(DailyPET, by = "dayIndex") %>% full_join(PETquant, by = "dayIndex") %>%
                    add_column(month = (dmy("1-09-2000") + 0:364))


        p4 = ggplot(data = SimPet, aes(x = month, group = 0)) + geom_line(aes(y = Measured, color = "Measured")) + ylab("PET [mm day\u207B]") +
                                geom_line(aes(y = PET, color = "Calculated"), size = 1.5, alpha = 0.5) +
                                scale_x_date(date_labels = "%B", date_breaks = "2 month")
        #-----
        #PET density      
        p5 = ggplot(data = IMSRain %>% filter(!is.na(pen))) + geom_density(aes(pen, color = "Measured"), bw = 0.05) + xlab("PET [mm day\u207B]") +
                                geom_density(data = SynthRain, aes(PET, color = "Calculated"), size = 1.5, alpha = 0.5, bw = 0.05)
}
    ggarrange(p1, p2,p3,p4,p5)
}

waterYearDayToMonth = function(day) {
    monthIndex = day %% 30
}