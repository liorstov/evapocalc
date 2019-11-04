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

CalculateProbabilities = function(IMSRain) {
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
    PWET = WetProb$counts / length(unique(IMSRain$year));

    #create probabilities for every annual day index
    Prob.Series = bind_cols("PWETOrig" = PWET, "PWAWOrig" = PWAW, "PWADOrig" = PWAD, dayIndex = 1:365) %>% replace_na(list(PWAWOrig = 0, PWADOrig = 0));
    Prob.Series = Prob.Series %>% add_column(month = month(dmy("1-10-2000") + Prob.Series$dayIndex));
    #smoothing
    window1 <- 50
    window2 = 40  
    Prob.Series = Prob.Series %>%
        mutate(PWET = stats::filter(stats::filter(Prob.Series$PWETOrig, rep(1 / window1, window1), method = "con", sides = 2, circular = 1), rep(1 / window2, window2), method = "con", circular = 1, sides = 2),
               PWAW = stats::filter(stats::filter(Prob.Series$PWAWOrig, rep(1 / window1, window1), method = "con", sides = 2, circular = 1), rep(1 / window2, window2), method = "con", circular = 1, sides = 2),
               PWAD = stats::filter(stats::filter(Prob.Series$PWADOrig, rep(1 / window1, window1), method = "con", sides = 2, circular = 1), rep(1 / window2, window2), method = "con", circular = 1, sides = 2))

   
    return(Prob.Series)
}
GetImsRain = function(station = 347700) {
    #stationElat = 347700;
    #stationSedom = 337000;
    con = odbcConnect("RainEvap")
    IMSRain = tbl_df(sqlQuery(con, paste("SELECT *,  year(time) as year,month(time) as month FROM data_dream.precip_daily where ((precip_daily.idstation = 347700))")));
    RODBC::odbcCloseAll();
    IMSRain = IMSRain %>% mutate(waterYear = ifelse(month %in% 10:12, year + 1, year)) %>%
                                 mutate(dayIndex = as.numeric(difftime(time, make_date(waterYear-1,10)))+1)

    IMSRain = IMSRain %>% filter(idstation == station)

    IMSRain = IMSRain %>% filter(year >= 1948 ,year <= 2017, measure > 0.1)
    

    return(IMSRain)
}
GenerateSeries = function(numOfyears = 69000, withEvapo = FALSE, IMSRain) {

    #    plot(pd, cex.axis = 1.5, cex.lab = 1.5, main = "Daily rainfall depth [mm]")
    DaysProb = CalculateProbabilities(IMSRain);
    #rain generator ---

    DepthLimit = 0.1
    #create the synthetic rain series and pick random values for occurance and possible amount from weibull
    #weibull parameters are taken for francesco

    tic()
    #create matrix with weibull values for depth
    weibull = matrix(rweibull(365 * numOfyears, shape = 0.4851, scale = 1.3458) + DepthLimit, nrow = numOfyears, ncol = 365);
    
    #matrix with random values 
    randMat = matrix(runif(365 * numOfyears), nrow = numOfyears, ncol = 365);
    SynthRain = matrix(0, nrow = numOfyears, ncol = 365, dimnames = list(1:numOfyears, 1:365))

    
    for (days in 2:ncol(SynthRain)) {
        SynthRain[, days] = as.numeric((SynthRain[, days - 1] & (randMat[, days] < DaysProb$PWAW[days])) | (!SynthRain[, days - 1] & (randMat[, days] < DaysProb$PWAD[days])));
    }
 
    SynthRain[which(SynthRain == 1)] = weibull[which(SynthRain == 1)];
    toc()

   
    #clear(c("weibull", "randMat"));


    if (withEvapo) {
        SynthRain = PETGen(SynthRain, IMSRain);
    }

    return(SynthRain) 
    #---
}

plotResults = function(SynthRain, IMSRain, DaysProb) {
    #graphic---
    SimRain = as_tibble(melt(SynthRain, value.name = "depth", varnames = c("year", "day")))
    numOfyears = length(unique(SimRain$year));

    #dividing the sim series to grups 
    SimRain$SeriesNumber = SimRain$year %/% length(unique(IMSRain$year));

    #wet day prob----
    SimRainDay = as_tibble((SynthRain > 0.1) * 1);
    SimRainDay = SimRainDay %>% add_column(series = (1:nrow(SimRainDay) - 1) %/% 69);
    SimRainDay = SimRainDay %>% group_by(series) %>% drop %>% dplyr::summarise_all(sum) / 69;
    bla1 = as_tibble(melt(SimRainDay[, 2:ncol(SimRainDay)])) %>% group_by(day = variable) %>%
                                dplyr::summarise(min = quantile(value, 0.05), median = quantile(value, 0.5), max = quantile(value, 0.95))
    quantile(SimRainDay[, 40])
    bla1= as_tibble(apply(SimRainDay, 2, FUN = quantile))
    bla2 = bla1 %>% gather(key = var_name, value = value) %>% spread(key = names(bla1), value = "value");

    SimRainDay = colSums(SimRainDay)
    SimRainDay = SimRain %>% filter(depth > 0.1) %>% group_by(SeriesNumber) %>% group_by(day) %>% arrange(SeriesNumber)
    SimRainDay = SimRain %>% filter(depth > 0.1) %>% group_by(day) %>% dplyr::summarise(n = n() / numOfyears, qu = quantiles()) 
    IMSday = IMSRain %>% group_by(day = dayIndex) %>% dplyr::summarise(n = n() / 69) 

    wetDayProb = SimRainDay %>% left_join(IMSday, by = "day", suffix = c(".sim", ".measure"))
    wetDayProb$n.measure = wetDayProb %>% dplyr::pull(n.measure) %>% replace_na(replace = 0)
    p1 =  ggplot(data = wetDayProb, aes(x = day, ymin = 0)) + geom_line(aes(y = n.measure, color = "Measured")) +
                geom_line(aes(y = n.sim, color = "Calculated"), size = 1.5) + ylab("Wet Probability [-]") + ggtitle(paste(numOfyears, "years"))
    #----
     
    #annual rain----
    simRainAnn = SimRain %>% filter(depth > 0.1) %>% group_by(year) %>% dplyr::summarise(annualRain = sum(depth), WetDays = n() )
    IMSRainAnn = IMSRain %>% filter(measure > 0.1) %>% group_by(waterYear) %>% dplyr::summarise(annualRain = sum(measure), WetDays = n())

    p2 = ggplot2::ggplot(data = simRainAnn, aes(annualRain)) + stat_density(aes(color = "Calculated"), size = 1.5, geom = "line", bw = 3) +
         stat_density(data = IMSRainAnn, aes(color = "Measured"), geom = "line", bw = 3) + xlim(0, 100) + xlab("Annual rain mm/year")
    mean(simRainAnn$WetDays)
    mean(IMSRainAnn$WetDays)
    #---
    #annual wet days
    p3 = ggplot2::ggplot(data = simRainAnn, aes(WetDays)) + stat_density(aes(color = "Calculated"), size = 1.5, geom = "line", bw = 1) +
         stat_density(data = IMSRainAnn, aes(WetDays, color = "Measured"), geom = "line", bw = 1) + scale_x_continuous(breaks = seq(0, 30, 5)) + xlab("# Annual wet days")
    #---
     print(ggarrange(p1,p2,p3)   )
}

waterYearDayToMonth = function(day) {
    monthIndex = day %% 30
}