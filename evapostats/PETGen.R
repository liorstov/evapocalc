require(lubridate)
require(moments)
require(plyr)

GetKforDay <- function(DayRow) {
    DayRow = as.numeric(DayRow);
    CurrentDay = DayRow[1];
    PrevDay = DayRow[3];
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
   
PETPerDay <- function(month, K) {
    dayCategoryIndex = which(K.month.table$K == K & K.month.table$month == month);
    meanPet = K.month.table$mean[dayCategoryIndex];
    std = K.month.table$std[dayCategoryIndex];
    rPet = K.month.table$rPen[dayCategoryIndex];
    SigmaCorrection = K.month.table$Corection[dayCategoryIndex];

    PET = meanPet - rPet * (PreviousDayGlobalVar - meanPet) + SigmaCorrection * std * (1 - rPet ^ 2) ^ 0.5;
    PreviousDayGlobalVar <<- PET;
    print(PreviousDayGlobalVar);
    return(PET);
}
PETGen <- function(RainSeries) {
    RainSeries = ReadFromGeoMateo("DB//EilotPenRain.csv");
    #Add column for previous day 
    RainSeries$prevDay = lag(RainSeries$rain, default = 0);

    #add column for previous day Pen
    RainSeries$PenLag = lag(RainSeries$Pen, default = 0);

    #add column for month
    RainSeries$month = lubridate::month(RainSeries$DateTime, label = FALSE)

    #k is the type of the day: 1:WAW,2:WAD,3:DAW,4DAD
    RainSeries$K = apply(RainSeries[, 2:4], 1, FUN = function(X) GetKforDay(X));

    #Prepare grouping by month and K
    RainSeries = RainSeries %>% group_by(K, month);

    #calculate statistics for each group
    K.month.table = RainSeries %>% dplyr::summarise(
                                             numElements = n(),
                                                 mean = mean(Pen),
                                                 std = sd(Pen),
                                                 skew = skewness(Pen),
                                                 rPen = cor(Pen, PenLag)                                                 
                                                )

    #remove NA std
    K.month.table$std[is.na(K.month.table$std)] = 0;

   #if number of measurments is small the set the corollation to 0.3
    K.month.table$rPen[which(K.month.table$numElements < 30)] = 0.3
    
    #calculate the gamma correction for the category
    K.month.table$Corection = apply(K.month.table[, c("skew", "rPen")], 1, FUN = function(X) GammaCorrection(X[1], X[2]));

    #adding PET to each day according to K value and month

    #Add column
    RainSeries$PET = NA;

    #The first day equal to the mean of its category
    firstDayIndex = which(K.month.table$K == RainSeries$K[1] & K.month.table$month == RainSeries$month[1]);
    RainSeries$PET[1] = K.month.table$mean[firstDayIndex];

    #calculate the rest of the days
    PreviousDayGlobalVar <<- RainSeries$PET[1];
    RainSeries$PET = apply(RainSeries[, c("month", "K", "prevDay")], 1, FUN = function(X) PETPerDay(X[1], X[2]));

    bla = RainSeries %>% dplyr::group_by(month) %>% dplyr::summarise(mean(PET))
}