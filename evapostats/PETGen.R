require(lubridate)
require(moments)
require(dplyr)
require(ggplot2)
require(reshape2)

GetKforDay <- function(CurrentDay, PrevDay) {
    if ((CurrentDay == 0) & (PrevDay == 0)) {
        return(1);
    }
    if (CurrentDay > 0 & PrevDay == 0) {
        return(2);
    }
    if (CurrentDay == 0 & PrevDay > 0) {
        return(3);
    }
    if (CurrentDay > 0 & PrevDay > 0) {
        return(4);
    }
}

ReadFromGeoMateo <- function(FileName) {
    RainSeries = na.omit(read.csv(FileName));
    RainSeries$DateTime = as.Date(RainSeries$DateTime, format = "%d/%m/%Y")
    return(RainSeries);
}

#Calculate Wilson Hilferty correction only ifskewness > 0.1 
GammaCorrection <- function(skew, rPen) {
    nt = rnorm(1);
    if (is.na(skew)) {
        return(nt);
    }
    else if (skew > 0.1) {
        gPet = (1 - rPen ^ 3) * skew / ((1 - rPen ^ 2) ^ 1.5);
        Correction = 2*((1 + gPet * nt / 6 - gPet ^ 2 / 36) ^ 3 - 1) / gPet;
        return(Correction);
    }
    else return(nt);
}
   
PETPerDay <- function(month, K, k.month.table) {
   
    dayCategoryIndex = head(which(k.month.table$K == K & k.month.table$month == month), 1);

    if (!length(dayCategoryIndex)) {
        return(PreviousDayGlobalVar);
    }
    meanPet = k.month.table$mean[dayCategoryIndex];
    std = k.month.table$std[dayCategoryIndex];
    rPet = k.month.table$rPen[dayCategoryIndex];
    SigmaCorrection = k.month.table$Corection[dayCategoryIndex];

    PET = meanPet - rPet * (PreviousDayGlobalVar - meanPet) + SigmaCorrection * std * (1 - rPet ^ 2) ^ 0.5;
    
    PreviousDayGlobalVar <<- PET;

    return(PET);
}


plotResults <- function() {


    #bind syntetic rain and PET
    bla = SyntRain %>% group_by(year) %>% dplyr::summarise(annualRain = sum(depth),
                                                            annualPET = sum(PET));

    bla$yeargroup = bla$year %/% 30;

    bla = bla %>% group_by(yeargroup) %>% dplyr::summarise(stdRain = sd(annualRain),
                                                            meanRain = mean(annualRain),
                                                            meanPET = mean(annualPET),
                                                            stdPET = mean(annualRain))
   

    blaIMS = IMSRain[which(IMSRain$year >= 1988),] %>% group_by(year) %>% dplyr::summarise(annualRain = sum(measure));
    IMSannualSD = sd(blaIMS$annualRain)
    IMSMeanAnnual = mean(blaIMS$annualRain)
    blaPET = RainSeries %>% group_by(year) %>% dplyr::summarise(annualPET = sum(measure));
    IMSPETstd = sd(blaPET$annualPET)
    IMSPETmean = mean(blaPET$annualPET)
    ##bind syntetic rain and PET
    #Synt = SyntRain %>% group_by(month);
    #Synt = Synt %>% dplyr::summarise(        
                                                 #meanPET = mean(PET),
                                                 #sumDepth = sum(depth),
                                                 #stdPET = sd(PET)
                                                #)
    #measured = RainSeries %>% group_by(month);
    #measured = measured %>% dplyr::summarise(meanPET = mean(Pen),
                                             #stdPET = sd(Pen))


    #stdHistogram rain
    ggplot2::ggplot(data = bla, aes(std)) + geom_histogram(aes(y = ..density..)) + geom_density(aes(color = "Blue"), show.legend = FALSE) +
        geom_vline(aes(xintercept = IMSannualSD, color = "red")) + geom_vline(aes(xintercept = density(bla$std)$x[which.max(density(bla$std)$y)], color = "Blue"), linetype = "dashed") +
        labs(title = "STD histogram for 30 yr chunks\nMeasured STD is 13.6" ) 
    mean(bla$std)

    theme_set(theme_gray(base_size = 18 , axis.text.x = element_text(size = 5)))
    geomTextSize =6;
    #mean histogram rain
    PickMAR = density(bla$meanRain)$x[which.max(density(bla$meanRain)$y)]
    ggplot2::ggplot(data = bla, aes(meanRain)) + geom_density(aes(color = "Calculated")) +
        geom_vline(aes(xintercept = IMSMeanAnnual, color = "Measured")) + geom_vline(aes(xintercept = PickMAR, color = "Calculated"), linetype = "dashed") +
        geom_text(aes(x = PickMAR + 0.3, label = paste(round(PickMAR,1), " [mm / yr]"), y = 0.04, color = "Calculated"), angle = 90, size = geomTextSize) +
        geom_text(aes(x = IMSMeanAnnual + 0.3, label = paste(round(IMSMeanAnnual,1), " [mm / yr]"), y = 0.04, color = "Measured"), angle = 90, size = geomTextSize) +
        labs(x = "Mean annual rainfall  [mm/yr]")

    PickSTD = density(bla$stdRain)$x[which.max(density(bla$stdRain)$y)]
    ggplot2::ggplot(data = bla, aes(stdRain)) + geom_density(aes(color = "Calculated")) +
        geom_vline(aes(xintercept = IMSannualSD, color = "Measured")) + geom_vline(aes(xintercept = density(bla$stdRain)$x[which.max(density(bla$stdRain)$y)], color = "Calculated"), linetype = "dashed") +
        geom_text(aes(x = PickSTD + 0.3, label = paste(round(PickSTD, 1), " [mm / yr]"), y = 0.04, color = "Calculated"), angle = 90, size = geomTextSize) +
        geom_text(aes(x = IMSannualSD + 0.3, label = paste(round(IMSannualSD, 1), " [mm / yr]"), y = 0.04, color = "Measured"), angle = 90, size = geomTextSize) +
        labs(x = "STD of Mean annual rainfall [mm/yr]")

    #petcomparison   
    ggplot2::ggplot(data = bla, aes(meanPET)) + geom_density(aes(color = "Calculated")) +
        geom_vline(aes(xintercept = IMSPETmean, color = "measured PET")) +
        labs(x = "Annual PET [mm/year]")
}

SynteticPetGen <- function(k.month.table, SyntRain) {

    SyntRain = as_tibble(SyntRain);

    ##rain as hydro year
    SyntRain$month = lubridate::month(lubridate::as_date(SyntRain$day, origin = dmy("01/01/1970") - 1))
    
    #Add column for previous day 
    SyntRain$prevDayRain = lag(SyntRain$depth, default = 0);

    #k is the type of the day: 1:DAD,2:WAD,3:DAW,4:WAW
    SyntRain$K = apply(SyntRain[, c("depth", "prevDayRain")], 1, FUN = function(X) GetKforDay(X[1], X[2]));

    #Add column
    SyntRain$PET = NA;

    #The first day equal to the mean of its category
    firstDayIndex = which(k.month.table$K == SyntRain$K[1] & k.month.table$month == SyntRain$month[1]);
    SyntRain$PET[1] = k.month.table$mean[firstDayIndex];

    #calculate the rest of the days
    PreviousDayGlobalVar <<- SyntRain$PET[1];
    SyntRain$PET = apply(SyntRain[, c("month", "K")], 1, FUN = function(X) PETPerDay(X[1], X[2], k.month.table));

    write.csv(SyntRain, "SyntRainPetEilat.csv")

    return(SyntRain)
}

PETGen <- function(SyntRain) {
    con = odbcConnect("RainEvap")
    RainSeries = tbl_df(sqlQuery(con, "SELECT *, measure as pen, year(time) as year,month(time) as month ,dayofyear(time) as dayIndex FROM data_dream.pet_daily where idstation = 347704 and isnull(measure) = 0"));

    RainSeries$rain = RainSeries %>% left_join(IMSRain, by = "time") %>% pull(measure.y) %>% replace_na(replace = 0)

    #Add column for previous day 

    ##becasue raw data is mean of 24 hours, still multiply by 24 is inaccurate
    #RainSeries$Pen = RainSeries$Pen * 24;

    RainSeries$prevDayRain = lag(RainSeries$rain, default = 0);

    #add column for previous day Pen
    RainSeries$PenLag = lag(RainSeries$pen, default = 0);

    ##add column for month
    #RainSeries$month = lubridate::month(RainSeries$DateTime, label = FALSE)

    #k is the type of the day: 1:DAD,2:WAD,3:DAW,4:WAW
    RainSeries$K = apply(RainSeries[, c("rain", "prevDayRain")], 1, FUN = function(X) GetKforDay(X[1], X[2]));

    #Prepare grouping by month and K
    RainSeries = RainSeries %>% group_by(K, month);

    #calculate statistics for each group
    #This table will supply the reference of penman evaporation for a specific K and month
    K.month.table = RainSeries %>% group_by(K, month) %>% dplyr::summarise(
                                                             numElements = n(),
                                                                 mean = mean(pen),
                                                                 std = sd(pen),
                                                                 skew = skewness(pen),
                                                                 rPen = cor(pen, PenLag)
                                                                )

    #remove NA std
    K.month.table$std[is.na(K.month.table$std)] = 0;

    #if number of measurments is small the set the corollation to 0.3
    K.month.table$rPen[which(K.month.table$numElements < 30)] = 0.3

    #calculate the gamma correction for the category
    K.month.table$Corection = apply(K.month.table[, c("skew", "rPen")], 1, FUN = function(X) GammaCorrection(X[1], X[2]));

    PET.rainfall.series = SynteticPetGen(K.month.table, SyntRain)
    return(PET.rainfall.series);
}