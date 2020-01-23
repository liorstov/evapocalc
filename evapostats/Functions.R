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

CalcGypsum <- function(raindata = SynthRain, duration, Depth = 100, thick = 2.5, wieltingPoint = 0.02, InitialCa = 0, initialSO4 = 0
                       , BulkDensity = 1.44, nArea = 1, FieldCapacity = 0.19, DustCa = 1.5, DustSO4 = 1.5, AETFactor = 0.6, RainFactor = 1, plotRes = TRUE, verbose = FALSE) {

    #if (!exists("b")) {
        #b = new(CSMCLASS);
    #}
    tic()
    list =  b$Calculate(raindata$rain*RainFactor, raindata$PET, duration, Depth, thick, wieltingPoint, InitialCa, initialSO4, BulkDensity, nArea, FieldCapacity,
                   DustCa, DustSO4, AETFactor, verbose);
    toc()
    list$thick = thick;
    #colnames(list$WD) = c("day", "rain", "AET", "WD", "soilMoisture", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20")
    #wd as Tibble
    DayWD = as_tibble(list$WD);
    list$Index30 = comp2Depth(which.min(abs(list$CompWash - (list$totalRain * 0.3))),thick = thick);
    list$Index03 = comp2Depth(which.min(abs(list$CompWash - (list$totalRain * 0.03))), thick = thick);
    list$WDp80 = round(quantile(DayWD$WD, 0.8),2);
    list$meanWD = round(mean(DayWD$WD), 2);

    #seasonal indices
    Seasonal = DayWD %>% mutate(year = day %/% 365) %>% group_by(year) %>% summarise(WD = max(WD))
    list$SMeanWD = round(mean(Seasonal$WD), 2);
    list$SWDp80 = round(quantile(Seasonal$WD, 0.8),2);
    if (plotRes) {
         plotSoilResults(list);
    }
    list$WD = DayWD %>% mutate(year = day %/% 365);
    list$waterBalance = list$totalRain - sum(list$AETLoss);
  return(list)  
}

comp2Depth = function(comp, thick) {
    return((comp-1+0.5)*thick)
}
plotSoilResults = function(res) {
    resLeach = tibble(leachate = res$CompWash, wetZone = res$WetZone / duration) %>% rowid_to_column %>%
                        mutate(WaterFluxOut = leachate / res$totalRain, depth = (rowid - 0.5) * res$thick, gypsum = res$gypsum)

    WP1 = ggplot(resLeach) + geom_bar(aes(x = depth, y = WaterFluxOut * 100), fill = "navyblue", stat = "identity", show.legend = FALSE) + scale_x_reverse(expand = c(0, 0.0)) + coord_flip() +
                    scale_y_continuous(name = "CaSO4 mEq/100g soil", position = "bottom", expand = c(0, 0.0), breaks = waiver()) +
                    theme(axis.text.x = element_text(size = 20, angle = 0, hjust = 1)) +
                    geom_line(aes(x = depth, y = gypsum, color = paste("gypsum =", resLeach$depth[which.max(res$gypsum)])), size = 1.3) +
                    geom_point(aes(x = depth, y = gypsum, color = paste("gypsum =", resLeach$depth[which.max(res$gypsum)])), size = 2) +
                                geom_vline(aes(xintercept = res$SMeanWD, color = paste("MeanSWD = ", res$SMeanWD))) +
                                geom_vline(aes(xintercept = res$Index03, color = paste("Index3 = ", res$Index03))) +
                                geom_vline(aes(xintercept = res$SWDp80, color = paste("SWDp80 = ", res$SWDp80)))


  
    #WP4 = ggplot() + annotate("text", x = 4, y = 25, label = paste("AET:", res$AET,"totalRain:", res$totalRain, sep = "\n")) + theme_bw() +
    #theme(panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank())

    return(ggarrange(WP1));
  # return(ggarrange(WP1, WP2, WP3,WP4,nrow = 2,ncol = 2))

}

plotMoisture = function(res) {

    resLeach = res$WD[nrow(res$WD), 6:20] %>% gather() %>% rowid_to_column %>%
                        mutate(depth = (rowid - 0.5) * res$thick, moisture = value / res$thick / 1.3 * 100, kik = c(11, 6.5, 4, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
    
    ggplot(resLeach) + geom_line(aes(x = depth, y = moisture),show.legend = FALSE) + scale_x_reverse() + coord_flip() +
                    scale_y_continuous(name = "moisture in soil ?g", position = "bottom", expand = c(0, 1), breaks  = seq(0,12,0.5)) +
                    theme(axis.text.x = element_text(size = 20, angle = 0, hjust = 1)) 

    return(rmsd(resLeach$moisture, resLeach$kik))
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