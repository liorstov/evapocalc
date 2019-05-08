require(lubridate)
require(moments)
require(plyr)

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
   
PETPerDay <- function(month, K) {
   
    dayCategoryIndex = head(which(K.month.table$K == K & K.month.table$month == month),1);

    if (!length(dayCategoryIndex)) {
        return(PreviousDayGlobalVar);
    }
    meanPet = K.month.table$mean[dayCategoryIndex];
    std = K.month.table$std[dayCategoryIndex];
    rPet = K.month.table$rPen[dayCategoryIndex];
    SigmaCorrection = K.month.table$Corection[dayCategoryIndex];

    PET = meanPet - rPet * (PreviousDayGlobalVar - meanPet) + SigmaCorrection * std * (1 - rPet ^ 2) ^ 0.5;
    
    PreviousDayGlobalVar <<- PET;

    return(PET);
}
PETGen <- function(fileName = "DB//EilotPenRain.csv") {
    RainSeries = ReadFromGeoMateo(FileName = fileName);
    #Add column for previous day 
    RainSeries$prevDay = lag(RainSeries$rain, default = 0);

    #add column for previous day Pen
    RainSeries$PenLag = lag(RainSeries$Pen, default = 0);

    #add column for month
    RainSeries$month = lubridate::month(RainSeries$DateTime, label = FALSE)

    #k is the type of the day: 1:WAW,2:WAD,3:DAW,4DAD
    RainSeries$K = apply(RainSeries[, 2:4], 1, FUN = function(X) GetKforDay(X[1],X[3]));

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


    measured = RainSeries %>% group_by(month);
    measured = measured %>% dplyr::summarise(
                                             numElements = n(),
                                                 meanMeasured = mean(Pen),
                                                 std = sd(Pen),
                                                 skew = skewness(Pen),                                               
                                                )
    Synt = raindata[1:365000,] %>% group_by(month);
    Synt = Synt %>% dplyr::summarise(        numElements = n(),
                                                 meanSynt = mean(PET)
                                               
                                                )
    PETBind = cbind2(measured, Synt);
    PETBind = melt(PETBind[,c(1,3,8)], id.vars = "month")
    ggplot(data = PETBind, aes(x = month, y = value, group = variable, color = variable)) + geom_line() +
    scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
     labs(x = "Mean PET[mm/d]", y = "Month", title = "PET in Eilat \ncalculated vs measured") ;
    return(RainSeries);
}