require(tidyverse)
require(MASS)
require(actuar)
library(RODBC)
library(tictoc)
require(fitdistrplus)
require(egg)
library(lubridate)
require(pracma)
require(cowplot)


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
        mutate(PWET = MovingAvarage(Prob.Series$PWETOrig,window1,window2),
               PWAW = MovingAvarage(Prob.Series$PWAWOrig, window1, window2),
               PWAD = MovingAvarage(Prob.Series$PWADOrig, window1, window2))

   
    return(Prob.Series)
}
MovingAvarage = function(x, win1, win2 = NULL) {
    ret = rollapply(x, win1, mean, na.rm = TRUE, fill = NA, partial = TRUE)
    if (!is.null(win2)) {
        ret = rollapply(ret, win2, mean, na.rm = TRUE, fill = NA, partial = TRUE)
    }

    #ret = na_ma(x, k = win1, weighting = "simple");
    #if (!is.null(win2)) {
        #ret = na_ma(x, k = win2, weighting = "simple");
    #}
    return (ret)
}

MovingAvarageSD = function(x, win1, win2 = NULL) {
    ret = rollapply(x, win1, sd, na.rm = TRUE, fill = NA, partial = TRUE)
    if (!is.null(win2)) {
        ret = rollapply(ret, win2, sd, na.rm = TRUE, fill = NA, partial = TRUE)
    }

    #ret = na_ma(x, k = win1, weighting = "simple");
    #if (!is.null(win2)) {
    #ret = na_ma(x, k = win2, weighting = "simple");
    #}
    return(ret)
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
GenerateRainSeries = function(NumOfSeries = 1000, IMSRain, AnuualRain = 0, WetDays = 0)
{
    DaysProb = CalculateProbabilities(IMSRain);    
    DepthLimit = 0.1
    NumOfSeries <<- NumOfSeries;
    #weibull parameters are taken for francesco
    station <<- getStationParam(IMSRain$station[1])
    stationName <<- station$stationName;
    shape = station$shape;
    scale = station$scale;
    factor = 1;
    numOfyears <<- measuredYears * NumOfSeries
    #create the synthetic rain series and pick random values for occurance and possible amount from weibull
    #rain generator ---
    print("Calculating rain:")
    tic()
    #if altered rain series
    if (AnuualRain != 0 & WetDays != 0) {
        if (stationName == "Eilat") {
            AltWB = AlternativeWeibullE(AnuualRain / WetDays);
            factor = 4.743 * WetDays ^ -0.711
        } else {
            AltWB = AlternativeWeibullS(AnuualRain / WetDays);
            factor = 5.057 * WetDays ^ -0.596

        }
        shape = AltWB$shape;
        scale = AltWB$scale;
    }



    #create matrix with weibull values for depth
    weibull = matrix(rweibull(365 * numOfyears, shape, scale) + DepthLimit, nrow = numOfyears, ncol = 365);

    #matrix with random values 
    randMat = matrix(runif(365 * numOfyears)*factor, nrow = numOfyears, ncol = 365);
    SynthRain = matrix(0, nrow = numOfyears, ncol = 365, dimnames = list(1:numOfyears, 1:365))
    
    for (days in 2:ncol(SynthRain)) {
        SynthRain[, days] = as.numeric((SynthRain[, days - 1] & (randMat[, days] < DaysProb$PWAW[days])) | (!SynthRain[, days - 1] & (randMat[, days] < DaysProb$PWAD[days])));
    }

    #addd amounts
    SynthRain[which(SynthRain == 1)] = weibull[which(SynthRain == 1)];
    toc()

    #as matrix

    SynthRain = as_tibble(melt(SynthRain, value.name = "rain", varnames = c("year", "dayIndex")))
    #clear(c("weibull", "randMat"));

   
    #dividing the sim series to grups 
    SynthRain = SynthRain %>% mutate(SeriesNumber = (year - 1) %/% measuredYears + 1)

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
    IMSRainDay = IMSRain %>% mutate(wet = (rain > 0.1) * 1) %>% group_by(dayIndex) %>% summarise(observed = sum(wet) / measuredYears)
    #hist for simulated

    SimRainDay = SynthRain %>% mutate(wet = (rain > 0.1) * 1) %>% group_by(dayIndex, SeriesNumber) %>% summarise(density = sum(wet) / measuredYears)
    SimRainDay = SimRainDay %>% group_by(dayIndex) %>%
                summarise(Simulated = mean(density), min = quantile(density, 0.05), median = quantile(density, 0.5), max = quantile(density, 0.95)) %>%
                    add_column(month = (dmy("1-09-2000") + 0:364)) %>% left_join(IMSRainDay, by = "dayIndex") %>%
                       left_join(rainProb, by="dayIndex");   

    p1 = ggplot(data = SimRainDay, aes(x = month, group = 1)) + geom_line(aes(y = observed, color = "Measured")) + ylab("Wet Prob. [-]") + xlab("") +
                                geom_line(aes(y = Simulated, color = "Simulated"), size = 1.5) +
                               # geom_line(aes(y = PWET, color = "WET")) +
                                ggtitle(paste(numOfyears, "years for", stationName)) +
                                geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.3) +
                                   scale_x_date(date_labels = "%B", expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0.0));

    rainProbabilities = ggplot(data = SimRainDay, aes(x = month, group = 1)) +                               
                                geom_line(aes(y = PWET, color = "WET")) +
                                geom_line(aes(y = PWAW, color = "WAW")) +
                                geom_line(aes(y = PWAD, color = "WAD")) +
                                   scale_x_date(date_labels = "%B", expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0.0))+ ylab("Wet Probability [-]") + xlab("");


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
    plotLimit = max(simRainAnn$WetDays);
    histBreaksSize = 1;
    #histogram for observed
    IMSAnnWetDays = hist(IMSRainAnn$WetDays, breaks = seq(0, plotLimit, histBreaksSize), plot = 0) %>%
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
                                     scale_y_continuous(expand = c(0, 0.0)) + scale_x_continuous(breaks = seq(0, plotLimit, 5), expand = c(0, 0));
    print("40%")

   
    #--
   
    p4=p5=ggplot()
    if (withEvapo) {

        #PET per day----
        plotLimit = 365;
        histBreaksSize = 1;
        #hist for observed
        IMSPet = IMSRain %>% filter(!is.na(pen)) %>% group_by(dayIndex) %>% summarise(observed = mean(pen))
        petYears = IMSRain %>% filter(!is.na(pen)) %>% distinct(waterYear) %>% nrow()
        SynthRain$PETSeries = SynthRain$year %/% petYears;

    
                   
        #aggregating   
        SimPETDay = SynthRain %>% group_by(dayIndex) %>% summarise(Simulated = mean(PET, na.rm = TRUE), min = quantile(PET, 0.05, na.rm = TRUE), median = quantile(PET, 0.5, na.rm = TRUE), max = quantile(PET, 0.95, na.rm = TRUE)) %>%
                left_join(IMSPet, by = "dayIndex") %>% left_join(PETProb, by = "dayIndex") %>% mutate(month = (dmy("1-09-2000") + dayIndex)) %>%
                mutate(K = replace(K, K == 1, "Wet"), K = replace(K, K == 2, "DAW"), K = replace(K, K == 3, "DAD"))
        p4 = ggplot(data = SimPETDay, aes(x = month, group = 1)) + geom_line(aes(y = Simulated, color = "Simulated"), size = 1.5) +
            geom_line(aes(y = observed, color = "Measured")) + ylab("PET [mm / day]") + xlab("") +
            geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.3) +
            scale_x_date(date_labels = "%B", expand = c(0, 0), label = "") + scale_y_continuous(expand = c(0, 0.0)) 
          #  geom_line(aes(group = K, y = smoothMean, color = as.factor(paste(K))))

        p7 = ggplot(data = SimPETDay, aes(x = month, group = 1))  +
            ylab("PET STD [mm / day]") + xlab("") +
            scale_x_date(date_labels = "%B", expand = c(0, 0), label = "") + scale_y_continuous(expand = c(0, 0.0)) + 
            geom_line(aes(group = K, y = smoothSTD, color = as.factor(paste(K))))
        print("60%")

        #--
        #PET density----
        plotLimit = ceiling(max(SynthRain$PET));
        histBreaksSize = 0.1;

        #histogram for observed
        DailyPET = IMSRain %>% filter(!is.na(pen)) %>% pull(pen) %>%
                                                        hist(breaks = seq(0, plotLimit, histBreaksSize), plot = 0) %>% Hist2tibble() %>% rename(DailyPET = density)

        DailyPETWet = IMSRain %>% filter(!is.na(pen)) %>% addKvalue() %>% filter(K == 1) %>% pull(pen) %>%
                                                        hist(breaks = seq(0, plotLimit, histBreaksSize), plot = 0) %>% Hist2tibble() %>% rename(DailyPETWet = density)

        DailyPETDAW = IMSRain %>% filter(!is.na(pen)) %>% addKvalue() %>% filter(K == 2) %>% pull(pen) %>%
                                                        hist(breaks = seq(0, plotLimit, histBreaksSize), plot = 0) %>% Hist2tibble() %>% rename(DailyPETDAW = density)

        DailyPETDAD = IMSRain %>% filter(!is.na(pen)) %>% addKvalue() %>% filter(K == 3) %>% pull(pen) %>%
                                                                hist(breaks = seq(0, plotLimit, histBreaksSize), plot = 0) %>% Hist2tibble() %>% rename(DailyPETDAD = density)
        
        #histogram for simulated
        SimPet = SynthRain %>% group_by(PETSeries) %>%
                       summarise(a = list(hist(PET, breaks = seq(0, plotLimit, histBreaksSize), plot = 0)), n = n()) %>% pull(a) %>%
                             map(~Hist2tibble(.x)) %>% reduce(bind_rows)

        #deivide to summer and winter
        #simPetWinter = SynthRain %>% filter(dayIndex <= 365 / 2) %>%
                       #summarise(a = list(hist(PET, breaks = seq(0, plotLimit, histBreaksSize), plot = 0)), n = n()) %>% pull(a) %>%
                             #map(~Hist2tibble(.x)) %>% reduce(bind_rows) %>% group_by(breaks) %>%
                        #summarise(simulatedWin = mean(density), minWin = quantile(density, 0.05), medianWin = quantile(density, 0.5), maxWin = quantile(density, 0.95))

        #simPetSummer = SynthRain %>% filter(dayIndex > 365 / 2) %>% 
                         #summarise(a = list(hist(PET, breaks = seq(0, plotLimit, histBreaksSize), plot = 0)), n = n()) %>% pull(a) %>%
                             #map(~Hist2tibble(.x)) %>% reduce(bind_rows) %>% group_by(breaks) %>%
                            #summarise(simulatedSum = mean(density), minSum = quantile(density, 0.05), medianSum = quantile(density, 0.5), maxSum = quantile(density, 0.95))

        #simPetDAW = SynthRain %>% filter(K == 2) %>% 
                         #summarise(a = list(hist(PET, breaks = seq(0, plotLimit, histBreaksSize), plot = 0)), n = n()) %>% pull(a) %>%
                             #map(~Hist2tibble(.x)) %>% reduce(bind_rows) %>% group_by(breaks) %>%
                            #summarise(simulatedDAW = mean(density), minDAW = quantile(density, 0.05), medianDAW = quantile(density, 0.5), maxDAW = quantile(density, 0.95))

        #rainOnly = SynthRain %>% filter(rain > 0) %>%
                         #summarise(a = list(hist(PET, breaks = seq(0, plotLimit, histBreaksSize), plot = 0)), n = n()) %>% pull(a) %>%
                             #map(~Hist2tibble(.x)) %>% reduce(bind_rows) %>% group_by(breaks) %>%
                            #summarise(simulatedRain = mean(density), minRain = quantile(density, 0.05), medianRain = quantile(density, 0.5), maxRain = quantile(density, 0.95))

        print("80%")

      
        #aggregating 
        SimPetDensity = SimPet %>% group_by(breaks) %>% summarise(simulated = mean(density), min = quantile(density, 0.05), median = quantile(density, 0.5), max = quantile(density, 0.95)) %>%
                left_join(DailyPET, by = "breaks") %>% left_join(DailyPETWet, by = "breaks") %>% left_join(DailyPETDAD, by = "breaks") %>% left_join(DailyPETDAW, by = "breaks") 

        p5 = ggplot(SimPetDensity, aes(x = breaks)) + geom_line(aes(y = simulated, color = "Simulated"), size = 1.5) + ylab("pdf") +
                                geom_line(aes(y = DailyPET, color = "Measured")) +
                                geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.3) +
                                xlab("PET [mm / day]") +
                                scale_y_continuous(expand = c(0, 0.0)) + scale_x_continuous(limits = c(0, plotLimit), expand = c(0, 0));

        PETProbabilities = ggplot(SimPetDensity, aes(x = breaks)) + geom_line(aes(y = simulated, color = "Simulated"), size = 1.5) + ylab("pdf") +
                                geom_line(aes(y = DailyPET, color = "Measured")) +
                                geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.3) +
                                geom_line(aes(y = DailyPETWet, color = "Measured wet")) +
                                geom_line(aes(y = DailyPETDAD, color = "Measured DAD")) +
                                geom_line(aes(y = DailyPETDAW, color = "Measured DAW")) +
                                xlab("PET [mm / day]") +
                                scale_y_continuous(expand = c(0, 0.0)) + scale_x_continuous(limits = c(0, plotLimit), expand = c(0, 0));

        #--
        
    }

  
    p6 = ggplot(tibble(values1 = rweibull(100000, 0.2, 4.5), values2 = rweibull(100000, station$shape, station$scale))) +
            stat_density(aes(values1, color = paste("Weibull\nScale =", 3.292875, "\nShape = ", 0.6640522)), size = 1, geom = "line") +
            stat_density(aes(values2, color = paste("Weibull\nScale =", station$shape, "\nShape = ", station$scale)), size = 1, geom = "line")+
               scale_y_continuous(expand = c(0, 0.0)) + scale_x_continuous(trans = 'log10', expand = c(0, 0), labels = scales::trans_format("log10", math_format(10 ^ .x))) + labs(x = "", y = "pdf");

    

    
    toc()

   # return(ggarrange(p1, rainProbabilities, p2, p3, p6, p4,  p7,p5, PETProbabilities))
    return(ggarrange(p1 + theme(legend.position = "none"), p2 + theme(legend.position = "none"), p3 + theme(legend.position = "none"), p4 + theme(legend.position = "none"), p5 + theme(legend.position = "none")))
   
  
  

}

waterYearDayToMonth = function(day) {
    monthIndex = day %% 30
}

AlternativeWeibullE = function(mean) {
    fun = function(x) { x * gamma(1 + 1 / (0.2 * log(x) + 0.4257)) - mean };
    scale = fzero(fun, c(1, 30))[[1]];
    shape = 0.2 * log(scale) + 0.4257;
    return(list(scale = scale, shape = shape))
}
AlternativeWeibullS = function(mean) {
    fun = function(x) { x * gamma(1 + 1 / (0.2 * log(x) + 0.5042)) - mean };
    scale = fzero(fun, c(0.75, 30))[[1]];
    shape = 0.2 * log(scale) + 0.5042;
    return(list(scale = scale, shape = shape))
}

GenerateSeries = function(station = 347700, stationEvap = 347704, NumOfSeries = 1000, AnuualRain = 0, WetDays = 0, PETfactor = 1) {
    if (PETfactor != 1) { 
        PETfactor = ifelse(station == 347700, PETfactor * 0.0005, PETfactor * 0.0004)
    }
    IMSRain = GetImsRain(station , stationEvap );
    rainSeriesResults = GenerateRainSeries(NumOfSeries = NumOfSeries, IMSRain = IMSRain, AnuualRain = AnuualRain, WetDays = WetDays);
    PETresults = PETGen(rainSeriesResults$SynthRain, IMSRain, 30);
    SynthRain = rainSeriesResults$SynthRain;
    SynthRain$PET = PETresults$SynthPET * PETfactor;
    SynthRain$K = PETresults$K;
    PETProb = PETresults$PETProb;
    rainProb = rainSeriesResults$DaysProb;
    SynthRain = SynthRain %>% arrange(year, dayIndex) %>% dplyr::select(year,rain ,PET);
    return(SynthRain)
}