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

CalcGypsum <- function(raindata = SedomData[, 3], years = 10000, Depth = 200, thick = 5, wieltingPoint = 0.04, InitialCa = 0, initialSO4 = 0
                       , BulkDensity = 1.44, nArea = 1, FieldCapacity = 0.4, DustCa = 1.5, DustSO4 = 1.5, AETFactor = 4.8, RainFactor = 1, observedArray = observedCom, Getresults = FALSE) {
    b = new(CSMCLASS);
    list = b$Calculate(raindata*RainFactor, years, Depth, thick, wieltingPoint, InitialCa, initialSO4, BulkDensity, nArea, FieldCapacity,
                   DustCa, DustSO4, AETFactor);
    observedArray = observedArray[!is.na(observedArray)]
    observedArray = cbind.fill(observedArray, list$gypsum, fill = NA)
    #print(observedArray)
    if (Getresults) {
        obsCalc <<- observedArray;
        print(observedArray)
    }

    rm(list);
    return(data.frame(difference(observedArray), years, Depth, thick, wieltingPoint, InitialCa, initialSO4, BulkDensity, nArea, FieldCapacity,
                   DustCa, DustSO4, RainFactor, AETFactor));
}

RawData2Compartments <- function(data, compartmentThick) {
    data$compartment = as.integer(data$depth_roof) %/% compartmentThick + 1

    #o.5 is the resolution
    ddply(data, .(compartment), numcolwise(function(x) sum(x * 0.5 / compartmentThick) ))
    
}