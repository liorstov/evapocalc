require(lubridate)
require(moments)
require(dplyr)
require(ggplot2)
require(reshape2)

#wet = 1, DAW = 2, DAD = 3
GetKforDay <- function(CurrentDay, PrevDay) {
    
    if (CurrentDay > 0) {
        return(1);
    }
    if (CurrentDay == 0 & PrevDay > 0) {
        return(2);
    }
    if ((CurrentDay == 0) & (PrevDay == 0)) {
        return(3);
    }    
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

#stochastic auto regressive model
PETPerDay <- function(monthodDay, Kday,random, K.month.table) {   
  
   print(Kday)
    ret = K.month.table %>% filter(month == monthodDay, K == Kday) %>% dplyr::select(mean, std);
    
    return(qnorm(random, ret$mean, ret$std));
}





#main function. parameter is a synthetic rain series and measured rain series
PETGen <- function(SynthRain, IMSRain) {
    tic()
    print("calc PET")   
    IMSPenOnly = IMSRain %>% filter(!is.na(pen))
    #Add column for previous day 
    IMSPenOnly$prevDayRain = lag(IMSPenOnly$rain, default = 0);

    #add column for previous day Pen
    IMSPenOnly$PenLag = lag(IMSPenOnly$pen, default = 0);    

    #k is the type of the day: 1:DAD,2:W,3:DAW,4:WAW
    IMSPenOnly$K = apply(IMSPenOnly[, c("rain", "prevDayRain")], 1, FUN = function(X) GetKforDay(X[1], X[2]));

    #Prepare grouping by day and K

    #calculate statistics for each group
    #This table will supply the reference of penman evaporation for a specific K and month
    K.month.table = IMSPenOnly %>% group_by(K, dayIndex) %>% dplyr::summarise(
                                                             numElements = n(),
                                                                 mean = mean(pen),
                                                                 std = sd(pen)                                                                
                                                                )
    #adding full days represenation
    K.month.table = tibble(K = rep(1:3, each = 365), dayIndex = rep(seq(365), 3)) %>% left_join(K.month.table, by = c("dayIndex", "K"))

    #filling NA days with selected mean values
    mean.per.day = K.month.table %>% group_by(dayIndex) %>% summarise(meanDay = min(mean, na.rm = 1))
    K.month.table$mean = apply(K.month.table, 1, FUN = function(X) { replace_na(X[4], mean.per.day$meanDay[X[2]]) })

    #remove NA std
    K.month.table$std[is.na(K.month.table$std)] = 0;

    windowSize = 30;
    K.month.table = K.month.table %>% group_by(K) %>% mutate(smoothMean = movavg(mean, windowSize, type = "s"), smoothSTD = movavg(std, windowSize, type = "s"))
   

    #calculate the gamma correction for the category
   # K.month.table$Corection = apply(K.month.table[, c("skew", "rPen")], 1, FUN = function(X) GammaCorrection(X[1], X[2]));
    PET.rainfall.series = SynteticPetGen(K.month.table, SynthRain)
    toc()
    return(list(SynthPET = PET.rainfall.series$PET, PETProb = K.month.table));
}

#set PET accoring to category
SynteticPetGen <- function(K.month.table, SynthRain) {
    
    ##rain as hydro year
   
    #Add column for previous day 
    SynthRain$prevDayRain = lag(SynthRain$rain, default = 0);

    #k is the type of the day: 1:DAD,2:WAD,3:DAW,4:WAW
    SynthRain$Wet = (SynthRain$rain > 0.1) *1
    SynthRain$daw = (!SynthRain$Wet & SynthRain$prevDayRain)*2 
    SynthRain$dad = (!SynthRain$Wet & !SynthRain$prevDayRain) * 3
    SynthRain = SynthRain %>% mutate(K = Wet + daw + dad) %>% dplyr::select(-Wet, - daw, - dad, - prevDayRain)

    #add random values
    SynthRain$rand = runif(nrow(SynthRain))

    SynthRain = SynthRain %>% left_join(K.month.table %>% dplyr::select(K, dayIndex, mean = smoothMean, std = smoothSTD), by = c("dayIndex", "K")) %>% mutate(PET = qnorm(rand, mean, std))%>%mutate(PET = ifelse(PET < 0 ,0,PET))
   

    return(SynthRain)
}