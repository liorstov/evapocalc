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


#add k values
#wet = 1, DAW = 2, DAD = 3
addKvalue <- function(IMSTable) {
    #Add column for previous day 
    IMSTable$prevDayRain = lag(IMSTable$rain, default = 0);

    #add column for previous day Pen
    IMSTable$PenLag = lag(IMSTable$pen, default = 0);

    #k is the type of the day: 1:DAD,2:W,3:DAW,4:WAW
    IMSTable$K = apply(IMSTable[, c("rain", "prevDayRain")], 1, FUN = function(X) GetKforDay(X[1], X[2]));
    IMSTable = IMSTable %>% dplyr::select(-PenLag);
    return(IMSTable);
}

#main function. parameter is a synthetic rain series and measured rain series
PETGen <- function(SynthRain, IMSRain, windowSize = 30, windowSize2 = NULL)
 {
    tic()
    print("calc PET")   

   #add k values
    IMSPenOnly = IMSRain %>% filter(!is.na(pen)) %>% arrange(time) %>% addKvalue();

    #representing every day
    #then pasting 2d moving avarage for every 30 days
    K.month.table = tibble(K = rep(1:3, each = 365), dayIndex = rep(seq(365), 3)) %>%
                    mutate(smoothMean = map2_dbl(K, dayIndex, ~ smoothMean(.x, .y, IMSPenOnly, windowSize)), smoothSTD = map2_dbl(K, dayIndex, ~ smoothStd(.x, .y, IMSPenOnly, windowSize))) %>%
                    mutate(smoothMean = replace_na(smoothMean, 0), smoothSTD = replace_na(smoothSTD,0))
   
    PET.rainfall.series = SynteticPetGen(K.month.table, SynthRain)
    toc() 
    return(list(SynthPET = PET.rainfall.series$PET,K = PET.rainfall.series$K, PETProb = K.month.table));
}

smoothMean = function(Kval,day, table, WS) {
    return(table %>% filter(K == Kval, dayIndex %in% cyclicSeq(day, WS)) %>% summarise(penT = mean(pen, na.rm = TRUE)) %>% pull(penT))
}
smoothStd = function(Kval, day, table,WS) {
    return(table %>% filter(K == Kval, dayIndex %in% cyclicSeq(day, WS)) %>% summarise(penSTD = sd(pen, na.rm = TRUE)) %>% pull(penSTD))
}
cyclicSeq = function(num, size) {
    size = size / 2;
    k_from = max(mod(num - size, 365), 1);
    k_to = max(mod(num + size, 365), 1);
    if (k_from < k_to) {
        ks = k_from:k_to
    }
    else
    { ks = c(k_from:365, 1:k_to) }
    return(ks)
    
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

    SynthRain = SynthRain %>% left_join(K.month.table %>% dplyr::select(K, dayIndex, mean = smoothMean, std = smoothSTD), by = c("dayIndex", "K")) %>%
                mutate(PET = qnorm(rand, mean, std)) %>% mutate(PET = ifelse(PET < 0, 0, PET)) %>% mutate(PET = replace_na(PET,5))
   

    return(SynthRain)
}