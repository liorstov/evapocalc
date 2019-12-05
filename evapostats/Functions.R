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



rmsd <- function(MatrixOC) {
    sqrt(sum((MatrixOC[, 1] - MatrixOC[, 2]) ^ 2))
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

CalcGypsum <- function(raindata = SynthRain, duration, Depth = 200, thick = 5, wieltingPoint = 0.04, InitialCa = 0, initialSO4 = 0
                       , BulkDensity = 1.44, nArea = 1, FieldCapacity = 0.4, DustCa = 1.5, DustSO4 = 1.5, AETFactor = 1, RainFactor = 1, Getresults = FALSE) {

    tic();
   # Rcpp::sourceCpp('C:/Users/liorst/source/repos/evapocalc/Calcyp/CSM.cpp', verbose = TRUE);
    raindata = raindata %>% arrange(year, dayIndex) %>% filter(year <= duration);
    b = new(CSMCLASS);
    b$Calculate(raindata$rain, raindata$PET, duration, Depth, thick, wieltingPoint, InitialCa, initialSO4, BulkDensity, nArea, FieldCapacity,
                   DustCa, DustSO4, AETFactor);
    list = b$GetResults();

    rm(b);
    toc();
    return(list$gypsum)
    #print("asdasd");
    #monthAET <<- (list$month);
    #WDVector <<- list$WD;
    #Gyps <<- list$gypsum;
    ##observedArray = observedArray[!is.na(observedArray)];
    ##observedArray = cbind.fill(observedArray, list$gypsum, fill = NA);
    ###print(observedArray)
    ##if (Getresults) {
        ##obsCalc <<- observedArray;
        ##print(observedArray)
    ##}

    #rm(list);
    #return(data.frame(years, Depth, thick, wieltingPoint, InitialCa, initialSO4, BulkDensity, nArea, FieldCapacity,
                   #DustCa, DustSO4, RainFactor, AETFactor));
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