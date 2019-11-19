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
   
    #dayCategoryIndex = head(which(K.month.table$K == K & K.month.table$month == month), 1);

    #if (!length(dayCategoryIndex)) {
        #return(PreviousDayGlobalVar);
    #}
  
    #rPet = K.month.table$rPen[dayCategoryIndex];
    #SigmaCorrection = K.month.table$Corection[dayCategoryIndex];

    #PET = meanPet - rPet * (PreviousDayGlobalVar - meanPet) + SigmaCorrection * std * (1 - rPet ^ 2) ^ 0.5;
    
    #PreviousDayGlobalVar <<- PET;
   print(Kday)
    ret = K.month.table %>% filter(month == monthodDay, K == Kday) %>% dplyr::select(mean, std);
    
    return(qnorm(random, ret$mean, ret$std));
}


#set PET accoring to category
SynteticPetGen <- function(K.month.table, SynthRain) {
    tic()
    ##rain as hydro year
    #SynthRain$month = lubridate::month(lubridate::as_date(SynthRain$day, origin = dmy("01/09/1970") - 1))
    
    #Add column for previous day 
    SynthRain$prevDayRain = lag(SynthRain$rain, default = 0);

    #k is the type of the day: 1:DAD,2:WAD,3:DAW,4:WAW
    SynthRain$Wet = (SynthRain$rain > 0.1) *1
    SynthRain$daw = (!SynthRain$Wet & SynthRain$prevDayRain)*2 
    SynthRain$dad = (!SynthRain$Wet & !SynthRain$prevDayRain)*3
    SynthRain = SynthRain %>% mutate(K = Wet + daw + dad) %>% dplyr::select(-Wet, - daw, - dad, -prevDayRain)
    #SynthRain$K = apply(SynthRain[, c("rain", "prevDayRain")], 1, FUN = function(X) GetKforDay(X[1], X[2]));

    #Add column
    #SynthRain$PET = NA;

    ##The first day equal to the mean of its category
    
    #firstDayIndex = which(K.month.table$K == SynthRain$K[1] & K.month.table$month == SynthRain$month[1]);
    #SynthRain$PET[1] = K.month.table$mean[firstDayIndex];

    ##calculate the rest of the days
    #PreviousDayGlobalVar <<- SynthRain$PET[1];
    #SynthRain$PET = apply(SynthRain[, c("month", "K")], 1, FUN = function(X) PETPerDay(X[1], X[2], K.month.table));
    SynthRain$rand = runif(nrow(SynthRain))
   # SynthRain1 = SynthRain[1:3,]  %>% rowwise() %>% mutate(PET = PETPerDay(month, K, rand, K.month.table))
    SynthRain = SynthRain %>% left_join(K.month.table %>% dplyr::select(K, dayIndex, mean = smoothMean, std = smoothSTD), by = c("dayIndex", "K")) %>%
                        mutate(PET = qnorm(rand, mean, std))
    #SynthRain = SynthRain %>% mutate(PET = replace_na(PET, qnorm(rand, mean(K.month.table$mean), mean(K.month.table$std))))
    paste("PET calc time: ",toc())

    return(SynthRain)
}

# for plotting, should be placed somware else
plotResults <- function() {


    #bind syntetic rain and PET
    bla = SynthRain %>% group_by(year) %>% dplyr::summarise(annualRain = sum(depth))
                                                            

    bla$yeargroup = bla$year %/% 30;

    bla = bla %>% group_by(yeargroup) %>% dplyr::summarise(stdRain = sd(annualRain),
                                                            meanRain = mean(annualRain))


    blaIMS = IMSRain %>% group_by(year) %>% dplyr::summarise(annualRain = sum(measure));
    IMSannualSD = sd(blaIMS$annualRain)
    IMSMeanAnnual = mean(blaIMS$annualRain)
    blaPET = RainSeries %>% group_by(year) %>% dplyr::summarise(annualPET = sum(measure));
    IMSPETstd = sd(blaPET$annualPET)
    IMSPETmean = mean(blaPET$annualPET)
   


    #stdHistogram rain
    ggplot2::ggplot(data = bla, aes(std)) + geom_histogram(aes(y = ..density..)) + geom_density(aes(color = "Blue"), show.legend = FALSE) +
        geom_vline(aes(xintercept = IMSannualSD, color = "red")) + geom_vline(aes(xintercept = density(bla$std)$x[which.max(density(bla$std)$y)], color = "Blue"), linetype = "dashed") +
        labs(title = "STD histogram for 30 yr chunks\nMeasured STD is 13.6")
    mean(bla$std)

    ggplot2::ggplot(data = bla, aes(y = meanRain)) + geom_boxplot() +
    ggplot(data = blaIMS, aes(y = annualRain)) + geom_boxplot()

    theme_set(theme_gray(base_size = 18, axis.text.x = element_text(size = 5)))
    geomTextSize = 6;
    #mean histogram rain
    PickMAR = density(bla$meanRain)$x[which.max(density(bla$meanRain)$y)]
    ggplot2::ggplot(data = bla, aes(meanRain)) + geom_density(aes(color = "Calculated")) +
        geom_vline(aes(xintercept = IMSMeanAnnual, color = "Measured")) + geom_vline(aes(xintercept = PickMAR, color = "Calculated"), linetype = "dashed") +
        geom_text(aes(x = PickMAR + 0, label = paste(round(PickMAR, 1), " [mm / yr]"), y = 0.04, color = "Calculated"), angle = 90, size = geomTextSize) +
        geom_text(aes(x = IMSMeanAnnual + 0, label = paste(round(IMSMeanAnnual, 1), " [mm / yr]"), y = 0.04, color = "Measured"), angle = 90, size = geomTextSize) +
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

#main function. parameter is a synthetic rain series and measured rain series
PETGen <- function(SynthRain, IMSRain) {

    IMSPenOnly = IMSRain %>% filter(!is.na(pen))
    #Add column for previous day 
    IMSPenOnly$prevDayRain = lag(IMSPenOnly$rain, default = 0);

    #add column for previous day Pen
    IMSPenOnly$PenLag = lag(IMSPenOnly$pen, default = 0);

    ##add column for month
    IMSPenOnly$month = lubridate::month(IMSPenOnly$time, label = FALSE)

    #k is the type of the day: 1:DAD,2:W,3:DAW,4:WAW
    IMSPenOnly$K = apply(IMSPenOnly[, c("rain", "prevDayRain")], 1, FUN = function(X) GetKforDay(X[1], X[2]));

    #Prepare grouping by day and K

    #calculate statistics for each group
    #This table will supply the reference of penman evaporation for a specific K and month
    K.month.table = IMSPenOnly %>% group_by(K, dayIndex) %>% dplyr::summarise(
                                                             numElements = n(),
                                                                 mean = mean(pen),
                                                                 std = sd(pen),
                                                                 rPen = cor(pen, PenLag)
                                                                )
    K.month.table = tibble(K = rep(1:3, each = 365), dayIndex = rep(seq(365), 3)) %>% left_join(K.month.table, by = c("dayIndex", "K"))
    k.max.mean = K.month.table %>% filter(K == 3) %>% dplyr::select(dayIndex, maxMean = mean)
    
    K.month.table$mean = apply(K.month.table, 1, FUN = function(X) { replace_na(X[4], k.max.mean$maxMean[X[2]]) })

    #remove NA std
    K.month.table$std[is.na(K.month.table$std)] = 0;

    K.month.table$smoothMean = 0;
    K.month.table$smoothSTD = 0;
    K.month.table$smoothMean[which(K.month.table$K == 1)] = K.month.table %>% filter(K == 1) %>% mutate(K == 1, d = movavg(mean, 30, type = "s")) %>% pull(d)
    K.month.table$smoothMean[which(K.month.table$K == 2)] = K.month.table %>% filter(K == 2) %>% mutate(K == 2, d = movavg(mean, 30, type = "s")) %>% pull(d)
    K.month.table$smoothMean[which(K.month.table$K == 3)] = K.month.table %>% filter(K == 3) %>% mutate(K == 3, d = movavg(mean, 30, type = "s")) %>% pull(d)

    K.month.table$smoothSTD[which(K.month.table$K == 1)] = K.month.table %>% filter(K == 1) %>% mutate(K == 1, d = movavg(std, 30, type = "s")) %>% pull(d)
    K.month.table$smoothSTD[which(K.month.table$K == 2)] = K.month.table %>% filter(K == 2) %>% mutate(K == 2, d = movavg(std, 30, type = "s")) %>% pull(d)
    K.month.table$smoothSTD[which(K.month.table$K == 3)] = K.month.table %>% filter(K == 3) %>% mutate(K == 3, d = movavg(std, 30, type = "s")) %>% pull(d)


    #calculate the gamma correction for the category
   # K.month.table$Corection = apply(K.month.table[, c("skew", "rPen")], 1, FUN = function(X) GammaCorrection(X[1], X[2]));
    PET.rainfall.series = SynteticPetGen(K.month.table, SynthRain)
    return(PET.rainfall.series$PET);
}

