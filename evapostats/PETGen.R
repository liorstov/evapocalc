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
PETGen <- function(fileName = "..//DB//EilotPenRain.csv") {
    RainSeries = ReadFromGeoMateo(FileName = fileName);
    #Add column for previous day 

    #becasue raw data is mean of 24 hours, still multiply by 24 is inaccurate
    RainSeries$Pen = RainSeries$Pen * 24;

    RainSeries$prevDay = lag(RainSeries$rain, default = 0);

    #add column for previous day Pen
    RainSeries$PenLag = lag(RainSeries$Pen, default = 0);

    #add column for month
    RainSeries$month = lubridate::month(RainSeries$DateTime, label = FALSE)

    #k is the type of the day: 1:DAD,2:WAD,3:DAW,4:WAW
    RainSeries$K = apply(RainSeries[, 2:4], 1, FUN = function(X) GetKforDay(X[1],X[3]));

    #Prepare grouping by month and K
    RainSeries = RainSeries %>% group_by(K, month);

    #calculate statistics for each group
    #This table will supply the reference of penman evaporation for a specific K and month
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

    return(K.month.table);
}

plotResults <- function() {


    #bind syntetic rain and PET
    bla = SyntRain %>% group_by(year) %>% dplyr::summarise(annualRain = sum(depth),
                                                            annualPET = sum(PET));

    bla$yeargroup = bla$year %/% 30;

    bla = bla %>% group_by(yeargroup) %>% dplyr::summarise(std = sd(annualRain),
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
    mean( bla$std)
    
    #mean histogram rain
    ggplot2::ggplot(data = bla, aes(mean)) + geom_histogram(aes(y = ..density..)) + geom_density(aes(color = "red"), show.legend = FALSE) +
        geom_vline(xintercept = IMSMeanAnnual, color = "blue") + geom_vline(xintercept = density(bla$mean)$x[which.max(density(bla$mean)$y)], color = "red", linetype = "dashed") +
        labs(title = "annual mean histogram for 30 yr chunks\nMeasured annual mean is 21.3 ")

    #petcomparison
    PETBind = cbind(measured = measured, synt = Synt);
    PETBind = melt(PETBind, id.vars = "month")
    ggplot(data = PETBind, aes(x = month, y = value, group = variable, color = variable)) + geom_line() +
    scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
     labs(x = "Mean PET[mm/d]", y = "Month", title = "PET in Eilat \ncalculated vs measured") ;
}

SynteticPetGen <- function(k.month.table, SyntRain) {

    SyntRain = as_tibble(SyntRain);

    #rain as hydro year
    SyntRain$month = lubridate::month(lubridate::as_date(SyntRain$day, origin = dmy("01/09/1970") - 1))

    #Add column for previous day 
    SyntRain$prevDay = lag(SyntRain$depth, default = 0);

    #k is the type of the day: 1:DAD,2:WAD,3:DAW,4:WAW
    SyntRain$K = apply(SyntRain[, c("depth", "prevDay")], 1, FUN = function(X) GetKforDay(X[1], X[2]));

    #Add column
    SyntRain$PET = NA;

    #The first day equal to the mean of its category
    firstDayIndex = which(K.month.table$K == SyntRain$K[1] & K.month.table$month == SyntRain$month[1]);
    SyntRain$PET[1] = K.month.table$mean[firstDayIndex];

    #calculate the rest of the days
    PreviousDayGlobalVar <<- SyntRain$PET[1];
    SyntRain$PET = apply(SyntRain[, c("month", "K")], 1, FUN = function(X) PETPerDay(X[1], X[2]));

    write.csv(SyntRain, "SyntRainPetEilat.csv")

    return(SyntRain[,c("month"  ,   "Depth", "PET")])
}

