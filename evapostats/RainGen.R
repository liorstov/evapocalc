require(tidyverse)
require(MASS)
require(actuar)
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
    #Prob.Series = Prob.Series %>% add_column(month = month(dmy("1-09-2000") + Prob.Series$dayIndex));
    #smoothing
    window1 <- 50
    window2 = 40  
    Prob.Series = Prob.Series %>%
        mutate(PWET = stats::filter(stats::filter(Prob.Series$PWETOrig, rep(1 / window1, window1), method = "con", sides = 2, circular = 1), rep(1 / window2, window2), method = "con", circular = 1, sides = 2),
               PWAW = stats::filter(stats::filter(Prob.Series$PWAWOrig, rep(1 / window1, window1), method = "con", sides = 2, circular = 1), rep(1 / window2, window2), method = "con", circular = 1, sides = 2),
               PWAD = stats::filter(stats::filter(Prob.Series$PWADOrig, rep(1 / window1, window1), method = "con", sides = 2, circular = 1), rep(1 / window2, window2), method = "con", circular = 1, sides = 2))

   
    return(Prob.Series)
}
GetImsRain = function(station = 347700 ,stationEvap = 347704) {
    #stationElat = 347700;
    #stationElatEvap = 347704;
    #stationSedom = 337000;
    con = odbcConnect("RainEvap")
    IMSRain = tbl_df(sqlQuery(con, paste("SELECT *, idstation as station, year(time) as year,month(time) as month, measure as rain FROM data_dream.precip_daily where ((precip_daily.idstation = ",station,"))")));
    EvapSeries = tbl_df(sqlQuery(con, paste("SELECT *, measure as pen, year(time) as year,month(time) as month ,dayofyear(time) as dayIndex FROM data_dream.pet_daily where idstation = ", stationEvap, " and isnull(measure) = 0")));
    RODBC::odbcCloseAll();

    #expend to get full year representation
    IMSRain = IMSRain %>% full_join(tibble(time = seq(IMSRain$time[1], IMSRain$time[1]+365, by = '1 day')), by = "time") %>% arrange(time)

    if (nrow(EvapSeries)) {
        IMSRain = IMSRain %>% full_join(EvapSeries, by = "time") %>% dplyr::select(time, rain, pen, station);
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
    station = getStationParam(IMSRain$station[1])
    stationName <<- station$stationName;
    
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

    #dividing the sim series to grups 
    SynthRain$SeriesNumber = SynthRain$SeriesNumber = rep(1:NumOfSeries, each = measuredYears)

    return(list(SynthRain = SynthRain, DaysProb = DaysProb))
    #---
}

Hist2tibble = function(histogram) {
    histogram = histogram %>% modify_at(1, ~ tail(.x, -1))
    return (tibble(breaks = histogram$breaks, counts = histogram$counts, density = histogram$density))
}

plotResults = function(SynthRain, IMSRain, rainProb, PETProb, withEvapo = FALSE) {
    #graphic--- 
    print("Calculating statistics:",tic())

    #wet day prob----   
    plotLimit = 365;
    histBreaksSize = 1;
    #hist for observed
    IMSRainDay = IMSRain %>% mutate(wet = (rain > 0) * 1) %>% group_by(dayIndex) %>% summarise(observed = sum(wet) / measuredYears)
    #hist for simulated

    SimRainDay = SynthRain %>% mutate(wet = (rain > 0) * 1) %>% group_by(dayIndex, SeriesNumber) %>% summarise(density = sum(wet) / measuredYears)
    SimRainDay = SimRainDay %>% group_by(dayIndex) %>%
                summarise(Simulated = mean(density), min = quantile(density, 0.05), median = quantile(density, 0.5), max = quantile(density, 0.95)) %>%
                    add_column(month = (dmy("1-09-2000") + 0:364)) %>% left_join(IMSRainDay, by = "dayIndex") %>%
                       left_join(rainProb, by="dayIndex");   

    p1 = ggplot(data = SimRainDay, aes(x = month, group = 1)) + geom_line(aes(y = observed, color = "Measured")) + ylab("Wet Probability [-]") +
                                geom_line(aes(y = Simulated, color = "Simulated"), size = 1.5) +
                                geom_line(aes(y = PWET, color = "PWET"), size = 1) +
                                geom_line(aes(y = PWAW, color = "PWAW"), size = 1) +
                                geom_line(aes(y = PWAD, color = "PWAD"), size = 1) +
                                ggtitle(paste(numOfyears, "years for", stationName)) +
                                geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.3) +
                                   scale_x_date(date_labels = "%B", expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0.0));

    #rm(SimRainDay, IMSRainDay)
    print("20%")

    #--
     
    #annual rain----

    #histogram for observed
    plotLimit = 120;
    histBreaksSize = 5;

    #histogram for observed
    IMSRainAnn = IMSRain %>% filter(rain > 0.1) %>% group_by(waterYear) %>% summarise(annual = sum(rain), WetDays = n()) %>% filter(annual <= plotLimit)
    IMShistannual = hist(IMSRainAnn$annual, breaks = seq(0, plotLimit, histBreaksSize), plot = 0) %>%
                        Hist2tibble() %>% mutate(obsDensity = counts / measuredYears)
    #histogram for simulated
    simRainAnn = SynthRain %>% filter(rain > 0.1) %>% group_by(year, SeriesNumber) %>% summarise(annual = sum(rain), WetDays = n()) %>% filter(annual <= plotLimit)
    rainHist = simRainAnn %>% group_by(SeriesNumber) %>%
                       summarise(a = list(hist(annual, breaks = seq(0, plotLimit, histBreaksSize), plot = 0)), n = n()) %>% pull(a) %>%
                        map(~Hist2tibble(.x)) %>% reduce(bind_rows)

    #aggregating 
    rainHist = rainHist %>% group_by(breaks) %>% summarise(simulated = mean(density), min = quantile(density, 0.05), median = quantile(density, 0.5), max = quantile(density, 0.95)) %>%
                left_join(IMShistannual, by = "breaks")

    p2 = ggplot(rainHist, aes(x = breaks)) + geom_line(aes(y = simulated, color = "Simulated"), size = 1.5) + ylab("pdf") +
                                geom_line(data = IMShistannual, aes(y = density, color = "Measured")) +
                                geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.3) + xlab("Annual rain mm/year") +
                                        scale_y_continuous(expand = c(0, 0.0)) + scale_x_continuous( expand = c(0,0));
    

   
    #--
    #annual wet days----
    plotLimit = 40;
    histBreaksSize = 1;
    #histogram for observed
    IMSAnnWetDays = hist(IMSRainAnn$WetDays, breaks = seq(-2, plotLimit, histBreaksSize), plot = 0) %>%
                        Hist2tibble() %>% rename(observed = density)

    #histogram for simulated
    SimAnnWetDays = simRainAnn %>% group_by(SeriesNumber) %>%
                       summarise(a = list(hist(WetDays, breaks = seq(-2, plotLimit, histBreaksSize), plot = 0)), n = n()) %>% pull(a) %>%
                        map(~Hist2tibble(.x)) %>% reduce(bind_rows)

    #aggregating 
    SimAnnWetDays = SimAnnWetDays %>% group_by(breaks) %>% summarise(simulated = mean(density), min = quantile(density, 0.05), median = quantile(density, 0.5), max = quantile(density, 0.95)) %>%
                left_join(IMSAnnWetDays, by = "breaks")

    p3 = ggplot(SimAnnWetDays, aes(x = breaks)) + geom_line(aes(y = simulated, color = "Simulated"), size = 1.5) + ylab("pdf") +
                                geom_line(aes(y = observed, color = "Measured")) +
                                geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.3) + xlab("# Annual wet days") +
                                     scale_y_continuous(expand = c(0, 0.0)) + scale_x_continuous(expand = c(0, 0));
    print("60%")

   
    #--
   
    p4=p5=ggplot()
    if (withEvapo) {

        #PET per day----
        plotLimit = 365;
        histBreaksSize = 1;
        #hist for observed
        IMSPet = IMSRain  %>% filter(!is.na(pen)) %>% group_by(dayIndex) %>% summarise(observed = mean(pen))
        #hist for simulated
        SimPETDay = SynthRain %>% group_by(dayIndex, SeriesNumber) %>% summarise(PET = mean(PET)) 
                   
        #aggregating   
        SimPETDay = SimPETDay %>% group_by(dayIndex) %>% summarise(Simulated = mean(PET), min = quantile(PET, 0.05), median = quantile(PET, 0.5), max = quantile(PET, 0.95)) %>%
                left_join(IMSPet, by = "dayIndex") %>% left_join(PETProb, by = "dayIndex") %>% mutate(month = (dmy("1-09-2000") + dayIndex)) %>%
                mutate(K = replace(K, K == 1, "Wet"), K = replace(K, K == 2, "DAW"), K = replace(K, K == 3, "DAD"))
        p4 = ggplot(data = SimPETDay, aes(x = month, group = 1)) + geom_line(aes(y = Simulated, color = "Simulated"), size = 1.5) +
            geom_line(aes(y = observed, color = "Measured total")) + ylab("PET [mm day\u207B]") +
            geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.3) +
            scale_x_date(date_labels = "%B", expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0.0)) +
            geom_line(aes(group = K, y = smoothMean, color = as.factor(paste("Expected",K))))

        #--
        #PET density----
        plotLimit = 20;
        histBreaksSize = 0.1;

        #histogram for observed
        DailyPET = IMSRain %>% filter(!is.na(pen)) %>%  pull(pen) %>%
                                                        hist(breaks = seq(0, plotLimit, histBreaksSize), plot = 0) %>% Hist2tibble() %>% rename(observed = density)

        #histogram for simulated
        SimPet = SynthRain %>% group_by(SeriesNumber) %>%
                       summarise(a = list(hist(PET, breaks = seq(0, plotLimit, histBreaksSize), plot = 0)), n = n()) %>% pull(a) %>%
                             map(~Hist2tibble(.x)) %>% reduce(bind_rows)
        print("80%")

      
        #aggregating 
        SimPetDensity = SimPet %>% group_by(breaks) %>% summarise(simulated = mean(density), min = quantile(density, 0.05), median = quantile(density, 0.5), max = quantile(density, 0.95)) %>%
                left_join(DailyPET, by = "breaks")

        p5=  ggplot(SimPetDensity, aes(x = breaks)) + geom_line(aes(y = simulated, color = "Simulated"), size = 1.5) + ylab("pdf") +
                                geom_line(aes(y = observed, color = "Measured")) +
                                geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.3) + xlab("PET [mm day\u207B]") +
                                     scale_y_continuous(expand = c(0, 0.0)) + scale_x_continuous(expand = c(0, 0));
        #--
        
    }
    ggarrange(p1, p2, p3, p4, p5)
    toc()
}

waterYearDayToMonth = function(day) {
    monthIndex = day %% 30
}