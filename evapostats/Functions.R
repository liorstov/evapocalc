CalcGypsum <- function(rainHolocene = SynthRainE, rainPleistocene = NULL, duration, Depth = 150, thick = 5, wieltingPoint = opt.WP,
                         nArea = 1, FieldCapacity = opt.FC, DustGyp = 0.005, AETFactor = opt.AETF,
                        verbose = FALSE, rainCa = 35.58, dustFlux = opt.dust, rainSO4 = opt.sulfate, getWD = F, random = T, distribution = F, rainstat = F,withRunoff = T,withFC = T,PETFactor=1) {


    regionEilat = !as.logical(rainHolocene %>% pull(year) %>% tail(1) %% 70)
    rainCa = rainSO4;
    # if pleistocene series suppied
    if (!is.null(rainPleistocene) & duration > 15000) {
        print("pleistocene profile")
        rainPleistocene$PET = rainPleistocene$PET*PETFactor
        #scrumble  both series
        if (random) {
            print("scramble rain series")
            tic()
            randomYears = sample(max(rainHolocene$year), 10000)
            rainHolocene = rainHolocene %>% filter(year %in% randomYears)
            randomYears = sample(max(rainPleistocene$year), duration - 10000)
            rainPleistocene = rainPleistocene %>% filter(year %in% randomYears)

            #combine series
            raindata = rainPleistocene %>% bind_rows(rainHolocene)
            raindata = raindata %>% mutate(origyear = year, year = rep(1:duration, each = 365))
            toc()
        }
        #no scrumble required
        else {
            raindata = rainPleistocene %>% filter(year <= (duration - 10000)) %>% bind_rows(rainHolocene %>% filter(year <= 10000)) %>% mutate(origyear = year, year = rep(1:duration, each = 365))
        }

        #apply runoff function
        if (withRunoff) {
            #identify region by length of series.           
            print(paste("is Eilat: ", regionEilat))

            #according to linear model in excel  
            if (regionEilat) {
                raindata = raindata %>% mutate(origRain = rain, rain = ageRainfuncEilat(year, rain));
            } else {
                raindata = raindata %>% mutate(origRain = rain, rain = ageRainfuncSedom(year, rain));
            }
            
        }
        else {
            raindata = raindata %>% mutate(origRain = rain)
        }

        rainStat = getSeriesMeans(rainPleistocene)
    } else {
        print("holocene profile")
        if (random) {
            #scrumble series
            randomYears = sample(max(rainHolocene$year), duration)
            raindata = rainHolocene %>% filter(year %in% randomYears)
            raindata = raindata %>% mutate(origyear = year, year = rep(1:duration, each = 365))

        }
        else {
            raindata = rainHolocene
        }

        #daysDuration = getDaysDuration(raindata);
        rainStat = getSeriesMeans(raindata)
    }


    print("run model")
    tic()
    list = cppModule$Calculate(raindata$rain, raindata$PET, duration, Depth, thick, wieltingPoint, nArea, FieldCapacity,
                   DustGyp, AETFactor, verbose, dustFlux / 10000 / 365, rainCa, rainSO4, withFC);
    toc()
    list$thick = thick;

    if (rainstat) {
        list$rainStat = rainStat
    } else {
        list$YearDust = NULL;
        list$YearGyp = NULL;
        list$YearMaxGyp = NULL;
        list$YearSulfate = NULL;
        list$YearCa = NULL;
        list$gypDepth = NULL
    }
    tic()
    list$max = mean(rainStat$max, na.rm = T)
    list$sd = mean(rainStat$sd, na.rm = T)
    list$rainDays = mean(rainStat$n, na.rm = T);
    list$AnnualRain = mean(rainStat$annual);
    list$q99 = rainStat %>% pull(max) %>% quantile(0.99) %>% as.numeric
    list$exc = rainStat %>% tally(max > 50) %>% as.numeric
    list$PET = mean(rainStat$PET, na.rm = T);
    if (distribution) {
        WeibullScaleShape = raindata %>% filter(rain > 0.1) %>% pull(rain) %>% fitdist("weibull", method = "mle")
        list$WBShape = WeibullScaleShape$estimate[1]
        list$WBScale = WeibullScaleShape$estimate[2]
    }
    #yearly = raindata %>% group_by(year) %>% summarise(yearlyRain = sum(rain), RainDays = length(which(rain > 0)), meanYear = mean(rain), RMSD = sd(rain));
    #list$YearlyRain = yearly %>% pull(yearlyRain);
    #list$RainDays = yearly %>% pull(RainDays);
    #list$meanYear = yearly %>% pull(meanYear);
    #list$RMSD = yearly %>% pull(RMSD);

    DayWD = dplyr::as_tibble(list$WD) %>% rowid_to_column(var = "day") %>% mutate(WD = value, year = day %/% 365) %>% filter(WD > 0)
    #list$Index30 = comp2Depth(which.min(abs(list$CompWash - (list$totalRain * 0.3))), thick = thick);
    #list$Index03 = comp2Depth(which.min(abs(list$CompWash - (list$totalRain * 0.03))), thick = thick);
    #list$WDp80 = round(quantile(DayWD$WD, 0.8), 2);
    #list$meanWD = round(mean(DayWD$WD), 2);

    ##seasonal indices
    Seasonal = DayWD %>% group_by(year) %>% summarise(maxWD = max(WD), meanWD = mean(WD))
    #list$WD = Seasonal
    list$SMeanWD = round(mean(Seasonal$maxWD), 2);
    list$SWDp80 = round(quantile(Seasonal$maxWD, 0.8), 2);

    #list$WD = DayWD %>% mutate(day = day + 1, year = (day + 1) %/% 365 + 1);
    list$maxWD = max(DayWD$WD);

    #DayGyp = dplyr::as_tibble(list$gypsumDay);
    #list$gypsumDay = DayGyp %>% mutate(day = day + 1, year = (day + 1) %/% 365 + 1);
    list$waterBalance = list$totalRain - sum(list$AETLoss);
    list$duration = duration;
    list$FC = FieldCapacity;
    list$WP = wieltingPoint;
    list$DF = dustFlux;
    list$DGyp = DustGyp;
    list$RCa = rainCa;
    list$RSO4 = rainSO4;
    list$AETF = AETFactor;

   toc()
    if (getWD) {
        list$WD = Seasonal
    }
    else {
        list$WD = NULL
    }
    #  list$WD = DayWD %>% filter(WD != 0) %>% mutate(year = day %/% 365) %>% group_by(year) %>% summarise(max = max(WD), sd = sd(max)) %>% mutate(cent = year %/% 100) %>% group_by(cent) %>% summarise(mean = mean(max), p80 = quantile(max, 0.8), sd = sd(max), var = var(max));
   # list$WD = Seasonal
    return(list)
}

#return not rmsd but sqrDiff netween obs abd res
RectanglingResults = function(res, obs) {

    ObsDepth = GetGypsicHorizonDepth(obs)
    ObsMaxGyp = GetGypsicHorizonConc(obs)
    ObsBottom = length(obs)

    obs = obs %>% pull(correctedMean)
    ObsGyp = sumGypsum(obs, 5,ObsBottom);
    #calculate target function
    rmsdTable = res %>% map_dfr(.f = ~(tibble(duration = .x$duration, AnnualRain = .x$AnnualRain, rainDays = .x$rainDays, sulfate = .x$RSO4, WP = .x$WP, FC = .x$FC, AETF = .x$AETF, dustFlux = .x$DF, PeakConc = max(.x$gypsum), PeakDepth = comp2Depth(which.max(.x$gypsum), .x$thick),
        total = sumGypsum(.x$gypsum, .x$thick, ObsBottom), SWDp80 = .x$SWDp80, obsTop = ObsDepth$top, obsBottom = ObsDepth$bottom, ObsGyp ,PET = .x$PET,AET = .x$AET))) %>% rowid_to_column()
    return(rmsdTable)

}

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

joinRMSD = function(diff, N) {
    return(sqrt(sum(diff) / sum(N)))
}


rmsd <- function(cal, obs) {

    #make sure same length
    depth = min(length(obs), length(cal))
    cal = cal[1:depth];
    obs = obs[1:depth];

    sqrt(sum((cal - obs) ^ 2) / length(cal))
}
MSD = function(observed, measured) {
    sum(abs(observed - measured)) / length(observed)
}



sqrDiff <- function(observed, measured) {
    limit = min(length(observed),length(measured))
    sum(((observed %>% head(limit)) - (measured %>% head(limit))) ^ 2)
}
LOOCV = function(CV) {
    resCV = tibble(profile = character(), val = numeric(), dust = numeric(), rain = numeric())
    for (item in unique(CV$profile)) {
        res = CV %>% filter(profile != item) %>% group_by(rain, dust) %>% summarise(RMSD = joinRMSD(mean, comps), meanNRMSD = mean(normRMSD))
        res = res %>% arrange(RMSD) %>% head(1)
        res = CV %>% filter(profile == item, dust == res$dust, rain == res$rain)
        resCV = resCV %>% add_row(profile = item, val = res$normRMSD, dust = res$dust, rain = res$rain)
    }

    return(resCV)
}



getSeriesMeans = function(series) {
    series %>% group_by(year) %>% summarise(n = sum(rain > 0.1), annual = sum(rain), max = max(rain), sd = sd(rain), PET = sum(PET), meanPET = mean(PET))
}
getDaysDuration = function(series) {
    temp = series %>% filter(rain > 0.1) %>% mutate(duration = ifelse(dayIndex - lag(dayIndex) == 1, 1, 0)) %>% group_by(year, grp = with(rle(duration), rep(seq_along(lengths), lengths)))  %>% summarise(duration = sum(duration) + 1) %>% group_by(year) %>% summarise(maxDuration = max(duration), meanDuration = mean(duration))
    return(temp)
}
CalcGypsumDouble <- function(rainHolocene = SynthRainE, rainPleistocene = NULL, duration, Depth = 100, thick = 5, wieltingPoint = opt.WP,
                         nArea = 1, FieldCapacity = opt.FC, DustGyp = 0.005, AETFactor = opt.AETF, plotRes = FALSE,
                        verbose = FALSE, rainCa = 35.58, dustFlux = opt.dust, rainSO4 = opt.sulfate, getWD = F, random = T) {

    list = CalcGypsum(rainHolocene, rainPleistocene, duration, Depth, thick, wieltingPoint,
                         nArea, FieldCapacity, DustGyp, AETFactor, plotRes,
                        verbose, rainCa, dustFlux, rainSO4, getWD = T, random = T);

    listOpt = CalcGypsum(list$WD, rainPleistocene, duration, Depth, thick, wieltingPoint = opt.WP,
                         nArea, FieldCapacity = opt.FC, DustGyp, AETFactor = opt.AETF, plotRes,
                        verbose, rainCa, dustFlux = opt.dust, rainSO4 = opt.sulfate, getWD = F, random = F);

    list$WD = NULL;
    listOpt$WD = NULL;
    list$ref =
    return(list(calc = list, ref = listOpt))

}

GetSensitivityFromResults = function(resultsTable) {

    #get values of the reference computation
    resultsTable$refPeakConc = resultsTable$PeakConc %>% lead(1)
    resultsTable$refPeakDepth = resultsTable$PeakDepth %>% lead(1)
    resultsTable$reftotal = resultsTable$total %>% lead(1)

    #erase ref computation
    resultsTable = resultsTable %>% filter(rowid %% 2 != 0 )

    #relative values from ref
    resultsTable = resultsTable %>% mutate(Max_Concentration = (PeakConc - refPeakConc) / refPeakConc,
                    Gypsum_depth = (PeakDepth - refPeakDepth) / refPeakDepth,
                       Total_concentration = (total - reftotal) / reftotal,WHC = FC-WP);

    #build tibble with values for each parameter and calculate the normelise factor
    senstest = bind_rows(
    WHCSens = resultsTable %>% filter(FC == opt.FC, AETF == opt.AETF, sulfate == opt.sulfate, dustFlux == opt.dust, WP %in% seq.WP, FC %in% seq.FC) %>% mutate(param = "WHC", change = WHC / opt.WHC, normal = opt.WHC / (WHC - opt.WHC)),
    WPSens = resultsTable %>% filter(FC == opt.FC, AETF == opt.AETF, sulfate == opt.sulfate, dustFlux == opt.dust, WP %in% seq.WP) %>% mutate(param = "\u03B8r", change = WP / opt.WP, normal = opt.WP / (WP - opt.WP)),
    FCSens = resultsTable %>% filter(WP == opt.WP, AETF == opt.AETF, sulfate == opt.sulfate, FC %in% seq.FC) %>% mutate(param = "FC", change = FC / opt.FC, normal = opt.FC / (FC - opt.FC)),
    AETFSens = resultsTable %>% filter(FC == opt.FC, WP == opt.WP, sulfate == opt.sulfate, dustFlux == opt.dust, AETF %in% seq.AETF) %>% mutate(param = "AET.F", change = AETF / opt.AETF, normal = opt.AETF / (AETF - opt.AETF)),
    sulfateSens = resultsTable %>% filter(FC == opt.FC, AETF == opt.AETF, WP == opt.WP, dustFlux == opt.dust, sulfate %in% seq.rainSeq) %>% mutate(param = "sulfate", change = sulfate / opt.sulfate, normal = opt.sulfate / (sulfate - opt.sulfate)),
    dustFSens = resultsTable %>% filter(FC == opt.FC, AETF == opt.AETF, sulfate == opt.sulfate, WP == opt.WP, dustFlux %in% seq.dustSeq) %>% mutate(param = "dustFlux", change = dustFlux / opt.dust, normal = opt.dust / (dustFlux - opt.dust))
    )

    #calculate with normalisation factor
    senstest = senstest %>% mutate(normal = abs(normal), Max_Concentration = Max_Concentration * normal, Gypsum_depth = Gypsum_depth * normal, Total_concentration = Total_concentration * normal) %>% drop_na()

    senstest$param = factor(senstest$param, levels = c("\u03B8r", "FC","WHC", "AET.F", "dustFlux", "sulfate"))
    #plot
    senstest %>% filter(change != 1) %>% dplyr::select(change, Max_Concentration, Gypsum_depth, Total_concentration, Parameter = param) %>% gather("target", "value", - Parameter, - change) %>%
    ggplot(aes(y = value, x = factor(change), fill = Parameter)) + facet_wrap(target ~ ., nrow = 1) + geom_boxplot(outlier.shape = NA) + stat_boxplot(geom = 'errorbar', linetype = 1, width = 0.5, position = position_dodge(0.75)) + geom_point(stat = "summary", fun = "mean", shape = 5, size = 2, position = position_dodge(0.75)) + geom_point(position = position_dodge(0.75)) + geom_hline(aes(yintercept = 0), linetype = "longdash") +
    scale_y_continuous(name = "Rel sensitivity [-]") + theme(strip.text = element_text(size = 30), legend.text = element_text(size = 30), legend.title = element_text(size = 35)) + labs(color = "Parameter", x = "change from optimal") + scale_fill_brewer(type = "qual", palette = 7) + scale_x_discrete(labels = function(x) paste0(as.numeric(x) * 100, "%"))

}

loadObsProfiles = function() {
    observedProfiles<<- read_csv("DB\\MeasuredDataProfiles.csv") 
    Zel11 <<- observedProfiles %>% filter(str_detect(SiteName, "11"))
    Zel12 <<- observedProfiles %>% filter(str_detect(SiteName, "12"))
    Zel13 <<- observedProfiles %>% filter(str_detect(SiteName, "13"))
    Zel1 <<- observedProfiles %>% filter(str_detect(SiteName, "ZEL1$"))
    Zel2 <<- observedProfiles %>% filter(str_detect(SiteName, "ZEL2$"))
    T1.9 <<- observedProfiles %>% filter(str_detect(SiteName, "T1-9"))
    T1.10 <<- observedProfiles %>% filter(str_detect(SiteName, "T1-10"))
    T1.1 <<- observedProfiles %>% filter(str_detect(SiteName, "T1-1$"))
    T2.1 <<- observedProfiles %>% filter(str_detect(SiteName, "T2-1$"))
    T2.10 <<- observedProfiles %>% filter(str_detect(SiteName, "T2-10"))
    Talus1 <<- observedProfiles %>% filter(str_detect(SiteName, "Talus"))
    SH1 <<- observedProfiles %>% filter(str_detect(SiteName, "Sheoheret 1"))
    SH3 <<- observedProfiles %>% filter(str_detect(SiteName, "Sheoheret 3"))
}
comp2Depth = function(comp, thick) {
    return((comp - 1 + 0.5) * thick)
}
plotSoilResults = function(res, observed, noInfo = F) {
    
    observedT = tibble(Observed = observed$correctedMean) %>% rowid_to_column
    resLeach = tibble(Simulated = res$gypsum) %>% rowid_to_column %>%
                        mutate(depth = (rowid - 0.5) * res$thick) %>% left_join(observedT) %>% mutate(inMeasured = is.na(Observed)) %>%
                        gather("factor", "gypsum", - depth, - rowid, - inMeasured)
    
    WP1 = ggplot(resLeach, aes(x = depth, y = gypsum, fill = factor, alpha = (inMeasured))) + geom_bar(stat = "identity", position = "dodge") + scale_x_reverse(expand = expansion(0, 2), name = "Depth [cm]", breaks = seq(0, 150, 30)) +
                    scale_y_continuous(name = "gypsum [mEq/100g soil]", expand = c(0, 0)) +
                    theme(axis.text.x = element_text(size = 20, angle = 0, hjust = 1), legend.title = element_blank()) +
                    guides(fill = guide_legend(reverse = TRUE), alpha = F) +
                    scale_fill_manual(values = c( "grey45","brown2")) + scale_alpha_discrete(range = c(1, 0.5)) +  coord_flip()
    if (!noInfo) {
       WP1 =  WP1 +
        ggtitle(paste("duration = ", res$duration, "; daily dust flux = ", round(res$DF * 1000, 3), "mg/cm3/day; dust SO4 = ", res$DF, "mg/g; Rain SO4 = ", res$RSO4, "mg/L; annual rain = \n", round(res$AnnualRain), "  rain days = ", round(res$rainDays), ";\ttotal simulated gypsum = ", round(sumGypsum(res$gypsum, 5, length(observed))), ";     total observed gypsum = ", round(sumGypsum(observed$correctedMean, 5, length(observed))))) +
        geom_vline(aes(xintercept = res$SWDp80, color = paste("SWDp80 = ", res$SWDp80)), linetype = 2, size = 2) + scale_color_manual(values = "orange")

    }
    #geom_line(aes(x = depth, y = Ca, color = "Ca")) +
    #geom_line(aes(x = depth, y = SO4, color = "SO4")) +
    #geom_vline(aes(xintercept = res$SMeanWD, color = paste("MeanSWD = ", res$SMeanWD))) +
    #geom_vline(aes(xintercept = res$Index03, color = paste("Index3 = ", res$Index03))) +
   
    ## geom_bar(aes(x = depth, y = WaterFluxOut * 100), fill = "navyblue", stat = "identity", show.legend = FALSE)
    #WP4 = ggplot() + annotate("text", x = 4, y = 25, label = paste("AET:", res$AET,"totalRain:", res$totalRain, sep = "\n")) + theme_bw() +
    #theme(panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank())

    return((WP1));
    # return(ggarrange(WP1, WP2, WP3,WP4,nrow = 2,ncol = 2))

}
plotSoilResultsAgg = function(res, obs, title = NULL) {
    SWDp80 = res[[1]]$SWDp80
    res = res %>% transpose %>% map_depth(2, ~ rowid_to_column(tibble(value = .x))) %>%
            pluck("gypsum") %>% melt(id.vars = "rowid", value.name = "value") %>% group_by(rowid) %>%
            dplyr::summarise(min = quantile(value, 0.05), calculated = mean(value), max = quantile(value, 0.95)) %>% mutate(depth = (rowid - 0.5) * res[[1]]$thick) %>%
            mutate(observed = obs, profile = title)

    res = res %>% dplyr::select(-rowid) %>% gather("factor", "gypsum", - depth, - min, - max, - profile)
    WP1 = ggplot(res, aes(x = depth, y = gypsum, fill = factor), palette = "jco") + geom_bar(stat = "identity", position = "dodge2") + scale_x_reverse(expand = c(0, 0.0)) + coord_flip(ylim = NULL) +
                    scale_y_continuous(name = "gypsum mEq/100g soil",expand = c(0, 2)) +
                    theme(axis.text.x = element_text(size = 20, angle = 0, hjust = 1)) +
                    geom_errorbar(aes(x = depth + 1, ymin = min, ymax = max)) +
                    geom_vline(aes(xintercept = SWDp80, color = paste("SWDp80"))) + theme(legend.title = element_blank()) + facet_wrap(profile ~ .) + theme(strip.text.x = element_text(size = 30))



    return(list(WP1, res));

}
plotSoilResultsRMSD = function(res, obs, CalibArray) {
    rmsdTable = res %>% transpose %>% map_depth(2, ~ rowid_to_column(tibble(value = .x))) %>%
    pluck("gypsum") %>% map_dbl(.f = ~rmsd(.x, obs)) %>% as_tibble() %>% add_column(CalibArray)
    WP <<- ggplot(rmsdTable, aes(x = CalibArray, y = value)) + geom_line() + geom_point() + ggtitle(paste("cor = ", cor(rmsdTable$CalibArray, rmsdTable$value)))
    print(WP)
    return(rmsdTable);

}

groupByParam = function(res) {
    groupTable = res %>% group_by(sulfate, dustFlux) %>% summarise(totalConc = sum(total), sumSqrDiff = sum(sqrDiff), min = quantile(sqrDiff, 0.05), max = quantile(sqrDiff, 0.95), sd = sd(sqrDiff), n = n(), comps = sum(comps), bias = mean(bias), PeakDepth = mean(PeakDepth), PeakConc = mean(PeakConc), minRMSD = joinRMSD(min, comps), RMSD = joinRMSD(sumSqrDiff, comps), maxRMSD = joinRMSD(max, comps)) %>% ungroup()
    return(groupTable)
}


#target function for overall diff
biasConc = function(res, obs, thick) {
    limit = min(length(obs), length(res))
    abs(sumGypsum(res, thick,limit) - sumGypsum(obs, thick,limit))
}



#gypsum concentration in undivided profile
sumGypsum = function(value, thick,bottom) {
    value = value[1:bottom];
    #calculate soil mass
    #value = value %>% pull(caso4)*0.4
    totalMass = length(value) * thick * 1.44;
    totalGypInProfile = sum(value/100 * thick*1.44)

    return(totalGypInProfile/totalMass*100)
}

#targetFunction for pick depth
biasDepth = function(res, obs) {
    abs(obs - which.max(res))
}

#targetFunction for pick concentration
PickConc = function(res, obs) {
    (max(obs) - max(res))
}


plotSoilResultsSurface = function(CV, required) {

    res = CV %>% filter(profile %in% c(required) & rain %in% seq(1, 20, by = 1) & dust %in% seq(0, 100, by = 1)) %>% group_by(rain, dust) %>% summarise(meanRMSD = joinRMSD(mean, comps), sum = sum(n), meanNRMSD = mean(normRMSD)) %>% mutate(profile = ifelse(length(required) > 1, "all", required))
    optimal = res %>% arrange(meanNRMSD) %>% filter(dust < 10) %>% head(1)
    WP = ggplot(res, aes(x = rain, y = dust, z = log10(meanRMSD), fill = log10(meanRMSD))) + scale_fill_gradientn(colours = rainbow(4)) + scale_x_continuous(expand = c(0, 0.0)) + scale_y_continuous(expand = c(0, 0.0)) +
        geom_raster(interpolate = TRUE) + geom_contour(color = "black", binwidth = 0.04) + labs(x = "Sulfate in rain water [mL/L]", y = "Dust Flux [g/m2/yr]", fill = "log(RMSD)") + facet_wrap(profile ~ .) + theme(strip.text.x = element_text(size = 30))
    # geom_point(data = optimal, size = 20, shape = 1, color = "blue", stroke = 3) + geom_text(data = optimal, size = 6, label = "Optimal point", color = "blue", hjust = 0, nudge_x = 1)
    return(list(WP, res))

    CV %>% filter(rain > 10, rain < 20, dust > 0, dust < 40, n > 10) %>% ggplot(aes((diffs), (PickDepth), color = rain, size = dust, shape = profile)) + labs(x = "bias [meq/100g soil]", y = "Depth peak diff [cm]", color = "sulfate mg/L", size = "dust flux mg/m2/yr") + geom_point() + scale_color_gradientn(colours = rainbow(5))
    CV %>% filter(rain > 10, rain < 20, dust > 0, dust < 40, n > 10) %>% ggplot(aes((diffs), (PickConc), color = rain, size = dust, shape = profile)) + labs(x = "bias [meq/100g soil]", y = " Concentration peak diff[meq / 100 g soil]", color = "sulfate mg/L", size = "dust flux mg/m2/yr") + geom_point() + scale_color_gradientn(colours = rainbow(5))
    CV %>% filter(rain > 10, rain < 25, dust > 0, dust < 40, n > 1) %>% ggplot(aes((meanRMSD), abs(diffs) / comps, color = abs(diffs) / comps / meanRMSD, size = dust, shape = profile)) + labs(x = "RMSD", y = "abs(bias)[meq/100g soil]", color = "bias/RMSD", size = "dust flux mg/m2/yr") + geom_point() + scale_color_gradientn(colours = rainbow(5))
    bla = CV %>% filter(rain > 0, rain < 25, dust > 0, dust < 50, n > 1) %>% group_by(rain, dust) %>% summarise(meanRMSD = joinRMSD(mean, comps), bias = abs(sum(diffs)) / sum(comps))

    require(ecr);
    test = CV %>% ungroup() %>% dplyr::select(RMSD, bias) %>% t() %>% which.nondominated()
    unloadNamespace("ecr")

    CV$pareto = FALSE;
    CV$pareto[test] = T;
    CV %>% ggplot(aes(RMSD, bias, color = factor(pareto))) + labs(x = "RMSD[meq/100g soil]", y = "abs(bias)[meq/100g soil]", color = "level", size = "dust flux mg/m2/yr") + geom_point(size = 3, show.legend = FALSE) + guides(color = guide_legend(override.aes = list(size = 8))) + coord_cartesian()
    CV %>% filter(pareto == 1) %>% dplyr::select(sulfate, dustFlux) %>% gather() %>% ggplot(aes(y = value, x = key)) + geom_boxplot() + labs(y = "dust flux/sulfate conc.", x = "param") + geom_point() + scale_y_continuous(breaks = scales::extended_breaks(10))

}

calculatePareto = function(parameterTable) {

    parameterTable = parameterTable %>% psel(low(biasConc) * low(biasDepth), top = nrow(parameterTable))
    parameterTable = parameterTable %>% mutate(distanceOrigin = abs(biasConc - biasDepth))

    print(parameterTable %>% ggplot(aes(biasConc, biasDepth, alpha = factor(.level < 3), color = .level)) + labs(x = "RMSD[meq/100g soil]", y = "abs(bias)[meq/100g soil]", size = "dust flux mg/m2/yr") + geom_point(size = 3) + scale_color_gradientn(colours = rainbow(5, rev = T)) + ggtitle(paste(parameterTable$profile[1], " ", parameterTable$duration[1], " yr")) + scale_alpha_discrete(guide = F))

    return(parameterTable)
}
plotSurfaceDiff = function(CV, required) {

    res = CV %>% filter(dust > 0, dust < 100, rain > 0) %>% group_by(profile) %>% nest() %>% mutate(interpulate = map(data, ~ interpulate(.x$rain, .x$dust, abs(.x$diffs) / .x$meanRMSD))) %>% unnest(interpulate) %>% dplyr::select(-data)
    res = res %>% bind_rows(tibble(profile = "all", x = res$x, y = res$y, z = res$z))
    WP = ggplot(res, aes(x = x, y = y, z = z, fill = z)) + scale_fill_gradientn(colours = rainbow(5)) + scale_x_continuous(expand = c(0, 0.0)) + scale_y_continuous(expand = c(0, 0.0)) +
        geom_raster(interpolate = TRUE) + geom_contour(color = "black") + labs(x = "Sulfate in rain water [mg/L]", y = "Dust Flux [g/m2/yr]", fill = "bias/RMSD[-]") + theme(strip.text.x = element_text(size = 30)) + facet_wrap(profile ~ ., 2)

    interpulate(bla$rain, bla$dust, bla$pareto) %>% ggplot(aes(x, y, fill = as.logical(z))) + geom_raster(interpolate = T, show.legend = FALSE) + coord_cartesian(ylim = c()) + labs(x = "Sulfate in rain water [mg/L]", y = "Dust Flux [g/m2/yr]", fill = "pareto level")
    ggplot(bla %>% filter(pareto == T), aes(x = rain, y = dust)) +
        geom_point() + labs(x = "Sulfate in rain water [mL/L]", y = "Dust Flux [g/m2/yr]", fill = "log(RMSD)") + theme(strip.text.x = element_text(size = 30))

    # geom_point(data = optimal, size = 20, shape = 1, color = "blue", stroke = 3) + geom_text(data = optimal, size = 6, label = "Optimal point", color = "blue", hjust = 0, nudge_x = 1)
    return(list(WP, res))
}

plotMoisture = function(res, obs) {



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

plotAnimation = function(res) {
    bla = res$gypsumDay[, c(1, 5:15)] %>% gather("comp", "gyp", - day, - WD) %>% mutate(depth = (as.numeric(comp) - 0.5) * 5, WD = ifelse(WD == 2.5, NA, WD))
    p = ggplot(bla) + geom_bar(aes(y = depth, x = gyp), stat = "identity") + scale_x_continuous(name = "gypsum meq/100g soil", position = "top") + geom_hline(aes(yintercept = WD, color = "Wetting\ndepth"), size = 3) +
                    scale_y_reverse(name = "depth [cm]") + scale_color_manual(values = "blue") +
                   theme(axis.text.x = element_text(size = 25, angle = 0, hjust = 0, margin = margin(b = 5))) + transition_time(day) + ggtitle('day: {as.integer(frame_time)}\t\t year: {as.integer(frame_time/365)}')
    animate(p, nframe = 400, renderer = av_renderer(), height = 800, width = 1200)
}

interpulate = function(x, y, z, n = 40) {
    mat = interp(x, y, z, duplicate = "median", ny = n, nx = n)
    d2 <- as_tibble(melt(mat$z, na.rm = T))
    names(d2) <- c("x", "y", "z")

    d2$x <- mat$x[d2$x]
    d2$y <- mat$y[d2$y]

    return(d2)
}
plotRate = function(res) {
    res = RateTest %>% map_df(.f = ~tibble(rain = .x$YearlyRain, RainDays = .x$RainDays, meanYear = .x$meanYear, RMSD = .x$RMSD, DFlux = .x$DF, RSO4 = .x$RSO4, gyp = .x$YearGyp, dustSulfate = .x$YearDust * 96.06, Rainsulfate = .x$YearSulfate * 96.06, gypsub = gyp, totSulfate = dustSulfate + Rainsulfate)) %>% drop_na() %>% filter(rain < 100 & rain > 0)
    ggplot(res, aes(dustSulfate * 1000, rain, color = gypsub, z = gypsub)) + geom_point(size = 1, alpha = 1) + scale_color_gradientn(colours = rainbow(5, rev = T)) + labs(x = "sulfate [mg/yr]", y = "yearly rain [mm]", color = "Gypsum \naccumulation rate \n [meq/100gr soil /yr]")


    bla = interp(res$ra, res$dustSulfate * 1000, res$gypsub, duplicate = "median", ny = 60, nx = 60)
    d2 <- as_tibble(melt(bla$z, na.rm = T))
    names(d2) <- c("rain", "totSulfate", "gypsub")

    d2$rain <- bla$x[d2$rain]
    d2$totSulfate <- bla$y[d2$totSulfate]
    ggplot(d2, aes(totSulfate, rain, fill = gypsub, z = gypsub, color = gypsub)) + geom_raster() + geom_contour(color = "black") + scale_fill_gradientn(colours = rainbow(5, rev = T)) + scale_color_gradientn(colours = rainbow(5, rev = T)) +
        labs(x = "sulfate from dust[mg/yr]", y = "annual standard deviation [mm]", fill = "Gypsum \naccumulation rate \n [meq/100gr soil /yr]") +
        geom_point(aes(x = 0.004, y = 30, shape = "T1.10 - 0.0025", fill = 3000), size = 5, color = "white", stroke = 3) +
        geom_point(aes(x = 0.01, y = 30, shape = "T1.9 - 0.0016"), color = "white", size = 5) +
        geom_point(aes(x = 0.018, y = 27, shape = "zel11 - 0.0018"), color = "white", size = 5) +
        geom_point(aes(x = 0.076, y = 30, shape = "T1.1 - 0.0145"), color = "white", size = 5) +
        geom_point(aes(x = 0.0099, y = 27, shape = "zel1 - 0.006"), color = "white", size = 5) +
        geom_point(aes(x = 0.0099, y = 30, shape = "T2.1 - 0.012"), color = "white", size = 5) +
        geom_point(aes(x = 0.004, y = 30, shape = "T1.10 - 0.0025", color = 0.0025), size = 3) +
        geom_point(aes(x = 0.01, y = 30, shape = "T1.9 - 0.0016", color = 0.0016), size = 3) +
        geom_point(aes(x = 0.018, y = 27, shape = "zel11 - 0.0018", color = 0.0018), size = 3) +
        geom_point(aes(x = 0.076, y = 30, shape = "T1.1 - 0.0145", color = 0.0145), size = 3) +
        geom_point(aes(x = 0.0099, y = 27, shape = "zel1 - 0.006", color = 0.006), size = 3) +
        geom_point(aes(x = 0.0099, y = 30, shape = "T2.1 - 0.012", color = 0.012), size = 3)

    d2 %>% group_by(rain) %>% summarise(mean = mean(gypsub)) %>% ggplot(aes(rain, mean)) + geom_bar(stat = "identity") + labs(y = " Gypsum \naccumulation rate \n [meq/100gr soil /yr]", x = "yearly rain [mm]")
    res %>% ggplot(aes(rain, gypsub)) + geom_point() + geom_point(data = res %>% filter(DFlux == 6 & RSO4 == 15), aes(color = "dustFlux = 6g/m2/yr \n rain = 15 mg/L")) + labs(y = "Gypsum \naccumulation rate \n [meq/100gr soil /yr]", x = "yearly rain [mm]")
}

##this function get a simulated rain series and PET to K table. returns the PET for every day in the
##Simulated series
#PETSeries <- function(raindata, K.Pet) {
    #raindata[, 4] = lag(raindata[, 3], default = 0);
    #raindata$K = apply(raindata[, 3:4], 1, FUN = function(X) GetKforDay(X[1], X[2]));
    #raindata$month = (raindata[, 2] %/% 31) + 1;

    ##The first day equal to the mean of its category
    #firstDayIndex = which(K.month.table$K == raindata$K[1] & K.month.table$month == raindata$month[1]);
    #raindata$PET[1] = K.month.table$mean[firstDayIndex];

    ##calculate the rest of the days
    #PreviousDayGlobalVar <<- raindata$PET[1];
    #raindata$PET[1:365000] = apply(raindata[1:365000, 5:6], 1, FUN = function(X) PETPerDay(X[2], X[1]));

#}


meanRMSD = function(results) {
    results = results ^ 2;
    return(sqrt(mean(results)));
}

ageRainfuncEilat = function(age, rain) {
  
    #get runoff ratio according to age and intensity
    #runoff = (0.0017 * log(age) - 0.0117) * (rain * factor)
    runoff = (5E-07*age + 0.0539) * rain

    #apply only if bigger then 10000
    runoff = runoff * (age > 10000)
    return(rain-runoff);
}
ageRainfuncSedom = function(age, rain) {

    #get runoff ratio according to age and intensity
    #runoff = (0.0017 * log(age) - 0.0117) * (rain * factor)
    runoff = (3E-07*age + 0.031) * rain

    #apply only if bigger then 10000
    runoff = runoff * (age > 10000)
    return(rain - runoff);
}

dailyRaintoIntensity = function(rain, isRegionEilat) {
    return(ifelse(isRegionEilat, 2.38 * rain, 1.76 * rain))
}

Theta.fc = function(silt, clay) {
    sand = 100 - silt - clay;
    rosetta = soilptf::soilinfo.rosetta((tibble(sand, clay, silt)))

    temp = ((rosetta$n - 1) / rosetta$n) ^ (1 - 2 * rosetta$n)
    Sfc = (1 + temp) ^ ((1 - rosetta$n) / rosetta$n)
    Theta.fc = rosetta$theta_r + Sfc * (rosetta$theta_s - rosetta$theta_r)
    return(Theta.fc)
}

Theta.r = function(silt, clay) {
    sand = 100 - silt - clay;
    rosetta = soilptf::soilinfo.rosetta((tibble(sand, clay, silt)))
   
    return(rosetta$theta_r)
}
addElementToResults = function() {
    results = results %>% modify_depth(2, ~ list_modify(.x, PET = ifelse(is.na(.x$PET), NA, .x$PET)))
    results = results %>% modify_depth(2, ~ list_modify(.x, duration = 0))
    bla = results$zel1 %>% map_depth(1, pluck("q99")) %>% map_df(tibble) %>% mutate(length("<list>"))
}


#ggplot <- function(...) ggplot2::ggplot(...) + scale_fill_brewer(type = "qual", palette = 2) + scale_color_brewer(type = "qual", palette = 2) 

GetResultsInRectangle = function(res, obs,sizeInPercent) {
    obsTotal = sumGypsum(obs %>% pull(caso4)*0.4, 5,length(obs));
    obsDepth = GetGypsicHorizonDepth(obs);

    print(c(obsTotal, obsDepth))

    sizeInPercentMin = 1-(sizeInPercent) / 100
    sizeInPercentMax = 1+(sizeInPercent) / 100
    return(res %>% mutate(inObsRange = ifelse(dplyr::between(PeakDepth, obsDepth * sizeInPercentMin, obsDepth * sizeInPercentMax) & dplyr::between(total, obsTotal * sizeInPercentMin, obsTotal * sizeInPercentMax), T, F)))
}

GetGypsicHorizonDepth = function(obs) {
   # obs %>% arrange(desc(top)) %>% group_by(caso4) %>% slice(1) %>% tail(1) %>% pull(bottom);
    return(obs %>% group_by(correctedMean) %>% summarise(top = min(top), bottom = max(bottom)) %>% tail(1));
}
GetGypsicHorizonConc = function(obs) {
    # obs %>% arrange(desc(top)) %>% group_by(caso4) %>% slice(1) %>% tail(1) %>% pull(bottom);
    return(obs %>% group_by(correctedMean) %>% summarise(bottom = mean(bottom)) %>% tail(1) %>% pull(correctedMean));
}

roundByNumber = function(value, factor) {
    value %/% factor * factor;
}

getPETfactor =  function(station = 347700, wantedPET) {
    return(ifelse(station == 347700, wantedPET * 5e-04, wantedPET *
        4e-04))
}