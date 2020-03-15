GetRandomValue <- function(X) {
    # if its a rain day
    if (interval == 0) {
        #set random value for rain amount 
        rand = runif(1, min = Amountcdf[1, 2]);
        X = Amountcdf[tail(which(Amountcdf[, 2] < rand), 1)];

        #get new interval
        rand = runif(1, min = Intervalcdf[1, 2]);
        interval <<- Intervalcdf[tail(which(Intervalcdf[, 2] < rand), 1)];

    }
    else {
        interval <<- interval - 1;
    }
}



rmsd <- function(observed,measured) {
    sqrt(sum((observed-measured) ^ 2)/length(observed))
}

difference <- function(MatrixOC) {
    depthofMaxValue = round(mean(which(max(MatrixOC[, 2], na.rm = TRUE) == MatrixOC[, 2]), na.rm = TRUE)) * 5;
    depthofMaxValueObserved = round(mean(which(max(MatrixOC[, 1], na.rm = TRUE) == MatrixOC[, 1]), na.rm = TRUE)) * 5;
    if (depthofMaxValue == 5) {
        depthofMaxValue = NA;
    }
    #print(depthofMaxValueObserved - depthofMaxValue)
    return((depthofMaxValueObserved - depthofMaxValue));
}

CalcGypsum <- function(raindata = SynthRain, duration, Depth = 100, thick = 5, wieltingPoint = 0.013,
                         nArea = 1, FieldCapacity = 0.1, DustGyp = 0.005, AETFactor = 1, plotRes = TRUE,
                        verbose = FALSE, dustFlux = 0.0152 / 365, rainCa = 35.58, rainSO4 = 20) {

       tic()
    list =  b$Calculate(raindata$rain, raindata$PET, duration, Depth, thick, wieltingPoint, nArea, FieldCapacity,
                   DustGyp, AETFactor, verbose, dustFlux, rainCa, rainSO4);
    toc()
    list$thick = thick;
    #colnames(list$WD) = c("day", "rain", "AET", "WD", "soilMoisture", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20")
    #wd as Tibble
    DayWD = as_tibble(list$WD);
    DayGyp = as_tibble(list$gypsumDay);
    list$Index30 = comp2Depth(which.min(abs(list$CompWash - (list$totalRain * 0.3))),thick = thick);
    list$Index03 = comp2Depth(which.min(abs(list$CompWash - (list$totalRain * 0.03))), thick = thick);
    list$WDp80 = round(quantile(DayWD$WD, 0.8),2);
    list$meanWD = round(mean(DayWD$WD), 2);

    #seasonal indices
    Seasonal = DayWD %>% mutate(year = day %/% 365) %>% group_by(year) %>% summarise(WD = max(WD))
    list$SMeanWD = round(mean(Seasonal$WD), 2);
    list$SWDp80 = round(quantile(Seasonal$WD, 0.8),2);
  
    list$WD = DayWD %>% mutate(day = day + 1, year = (day + 1) %/% 365 + 1);
    list$gypsumDay = DayGyp %>% mutate(day = day + 1, year = (day + 1) %/% 365 + 1);
    list$waterBalance = list$totalRain - sum(list$AETLoss);
    list$duration = duration;
    list$FC = FieldCapacity;
    list$WP = wieltingPoint;
    list$DF = dustFlux;
    list$DGyp = DustGyp;
    list$RCa = rainCa;
    list$RSO4 = rainSO4;
    list$maxWD = max(DayWD$WD);

    if (plotRes) {
            plotSoilResults(list);
    }
    return(list)  
}

comp2Depth = function(comp, thick) {
    return((comp-1+0.5)*thick)
}
plotSoilResults = function(res,observed = NULL) {
    resLeach = tibble(leachate = res$CompWash, wetZone = res$WetZone / res$duration) %>% rowid_to_column %>%
                        mutate(WaterFluxOut = leachate / res$totalRain, depth = (rowid - 0.5) * res$thick, calculated = res$gypsum, observed) %>%
                        dplyr::select(depth, calculated, observed) %>% gather("factor", "gypsum", - depth)

    WP1 = ggplot(resLeach, aes(x = depth, y = gypsum, fill = factor)) + scale_x_reverse(expand = c(0, 0.0)) + coord_flip(ylim = NULL) +
                    scale_y_continuous(name = "gypsum mEq/100g soil", position = "bottom", expand = c(0, 0.0)) +
                    theme(axis.text.x = element_text(size = 20, angle = 0, hjust = 1)) +
                   geom_bar(stat = "identity", position = "dodge2") +
                    #geom_line(aes(x = depth, y = Ca, color = "Ca")) +
                    #geom_line(aes(x = depth, y = SO4, color = "SO4")) +
                    ggtitle(paste("duration = ", duration, "; daily dust flux = ", round(res$DF*1000,3), "mg/cm3/day; dust Ca = ", res$DCa, "mg/g; dust SO4 = ", res$DSO4, "mg/g; Rain Ca = ", res$RCa, "mg/L; Rain SO4 = ", res$RSO4,"mg/L")) +
                    #geom_vline(aes(xintercept = res$SMeanWD, color = paste("MeanSWD = ", res$SMeanWD))) +
                    #geom_vline(aes(xintercept = res$Index03, color = paste("Index3 = ", res$Index03))) +
                    geom_vline(aes(xintercept = res$SWDp80, color = paste("SWDp80 = ", res$SWDp80))) 
             ## geom_bar(aes(x = depth, y = WaterFluxOut * 100), fill = "navyblue", stat = "identity", show.legend = FALSE)
    #WP4 = ggplot() + annotate("text", x = 4, y = 25, label = paste("AET:", res$AET,"totalRain:", res$totalRain, sep = "\n")) + theme_bw() +
    #theme(panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank())

    return(ggarrange(WP1));
  # return(ggarrange(WP1, WP2, WP3,WP4,nrow = 2,ncol = 2))

}
plotSoilResultsAgg = function(res) {

    res = res %>% dplyr::select(-rowid) %>% gather("factor", "gypsum", -depth,-min,-max)
    WP1 = ggplot(res, aes(x = depth, y = gypsum, fill = factor)) + scale_x_reverse(expand = c(0, 0.0)) + coord_flip(ylim = NULL) +
                    scale_y_continuous(name = "gypsum mEq/100g soil", position = "bottom") +
                    theme(axis.text.x = element_text(size = 20, angle = 0, hjust = 1)) +
                   geom_bar(stat = "identity", position = "dodge2") +
                    geom_errorbar(aes(x = depth+1, ymin = min, ymax = max) )
   

    return(ggarrange(WP1));
    # return(ggarrange(WP1, WP2, WP3,WP4,nrow = 2,ncol = 2))

}
plotMoisture = function(res,obs) {


    
    #48 hours late holo from kiki - 3, 2, 1.5, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
    #12 late holo 11, 6.5, 4, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
    #fig 7-1 48Hours late holo 6, 4, 1.5, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
    #c(3.2, 2, 1.4, 1, 1, 1, 1, 1)
    moistPercentage = res$WD[nrow(res$WD), 6:10] %>% gather() %>% pull(value) %>% sum() / 50 / 1.3 * 100
    #return(abs(1-moistPercentage))
    resLeach = res$WD[nrow(res$WD), 6:10] %>% gather() %>% rowid_to_column %>%
                        mutate(depth = (rowid - 0.5) * res$thick, simulated = value / res$thick / 1.3 * 100) %>%
                        mutate(observed = obs) %>% dplyr::select(depth, observed, simulated) %>% gather("source", "moist", - depth)
    
    ggplot(resLeach, aes(x = depth, y = moist, fill = source)) + geom_bar(stat = "identity", position = "dodge2") + scale_x_reverse() + coord_flip() +
                    scale_y_continuous(name = "moisture in soil [\u03B8g] \n Kc = 1.2; WP = 0.013 cm3/cm3; FC = 0.1 cm3/cm3", position = "bottom", expand = c(0, 0.008)) +
                    theme(axis.text.x = element_text(size = 20, angle = 0, hjust = 1)) 

    return(rmsd(resLeach$sim, resLeach$kik))
}

plotAnimation= function(res, obs) {
    bla = res$gypsumDay[, c(1, 6:15)] %>% gather("comp", "gyp", - day) %>% mutate(depth = (as.numeric(comp) - 0.5) * 5)
    p = ggplot(bla) + geom_bar(aes(x = depth, y = gyp), stat = "identity") + scale_x_reverse() + coord_flip() +
                    scale_y_continuous(name = "gypsum meq/100g soil", position = "bottom") +
                    theme(axis.text.x = element_text(size = 20, angle = 0, hjust = 1)) + transition_time(day) + ggtitle('day: {frame_time}')
    animate(p, nframe = 200, renderer = av_renderer())
}

#this function get a simulated rain series and PET to K table. returns the PET for every day in the
#Simulated series
PETSeries <- function(raindata, K.Pet) {
    raindata[, 4] = lag(raindata[, 3], default = 0);
    raindata$K = apply(raindata[, 3:4], 1, FUN = function(X) GetKforDay(X[1], X[2]));
    raindata$month = (raindata[, 2] %/% 31) + 1;

    #The first day equal to the mean of its category
    firstDayIndex = which(K.month.table$K == raindata$K[1] & K.month.table$month == raindata$month[1]);
    raindata$PET[1] = K.month.table$mean[firstDayIndex];

    #calculate the rest of the days
    PreviousDayGlobalVar <<- raindata$PET[1];
    raindata$PET[1:365000] = apply(raindata[1:365000, 5:6], 1, FUN = function(X) PETPerDay(X[2],X[1]));
    
}

RawData2Compartments <- function(data, compartmentThick) {
    data$compartment = as.integer(data$depth_roof) %/% compartmentThick + 1

    #o.5 is the resolution
    ddply(data, .(compartment), numcolwise(function(x) sum(x * 0.5 / compartmentThick) ))
    
}