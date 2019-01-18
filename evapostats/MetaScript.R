
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
    depthofMaxValue = min(which(max(MatrixOC[, 2]) == MatrixOC[, 2])) * 5;
    depthofMaxValueObserved = min(which(max(MatrixOC[, 1]) == MatrixOC[, 1])) * 5;
    return(abs(depthofMaxValueObserved - depthofMaxValue));
    }

CalcGypsum <- function(raindata = EilatData[, 1], years = 100, Depth = 100, thick = 5, wieltingPoint = 0.04, InitialCa = 10, initialSO4 = 10
                       , BulkDensity = 1.44, nArea = 1, FieldCapacity = 0.4, DustCa = 500, DustSO4 = 500, AETFactor = 1, observedArray = observedCom, Getresults = TRUE) {
    b = new(CSMCLASS);
    list = b$Calculate(raindata, years, Depth, thick, wieltingPoint, InitialCa, initialSO4, BulkDensity, nArea, FieldCapacity,
                   DustCa, DustSO4, AETFactor);
    #cat(wieltingPoint, "      ", list$gypsum, '\n')
    observedArray[, 2] = list$gypsum;

    if (Getresults) {
        obsCalc <<- observedArray;
    }
    rm(list);
    return(data.frame(difference(observedArray), years, Depth, thick, wieltingPoint, InitialCa, initialSO4, BulkDensity, nArea, FieldCapacity,
                   DustCa, DustSO4, AETFactor));
    }

    # temp = CalcGypsum(years = 10000, DustCa = 5, DustSO4 = 5, AETFactor = 100, wieltingPoint = 0.04, FieldCapacity = 0.4)

require(Rcpp)
require(ggplot2)
require(reshape2)

#set wd
setwd("C:/Users/Lior/master/evapocalc/");
EilatData = read.csv("DB/rainseriestemp.csv");
Observed = read.csv("DB/measured.CSV");






#create observed comartments
observedCom = matrix(nrow=20, ncol = 2);
observedCom[1,1] = 0.72;
observedCom[2:8,1] = 2.44;
observedCom[8:20, 1] = 4.2;
colnames(observedCom) = c("Observed", "Calculated");


Rcpp::sourceCpp('C:/Users/Lior/master/evapocalc/Calcyp/CSM.cpp', verbose = TRUE);
temp = CalcGypsum(years = 100);
wiltingPointArray = seq(0.001, 0.04, 0.001);
DustFluxArray = seq(0.5, 5, 0.1);
AETArray = seq(from = 100, to = 300, length = 45);

#running differrent values
results = sapply(DustFluxArray, FUN = function(X) CalcGypsum(years = 10000, DustCa = X, DustSO4 = X, AETFactor = 100));

#evapotranspiration
results = sapply(AETArray, FUN = function(X) CalcGypsum(years = 10000, DustCa = 2, DustSO4 = 2, AETFactor = X));

#convert to plotable
results = as.data.frame.array(t(results))
for (I in seq(1:ncol(results))) {
    results[, I] = unlist(results[I]);
}

#sesitivity test for flux
dustplot = ggplot(data = results) + geom_smooth(se = FALSE, mapping = aes(x = results$DustCa, y = results$difference.observedArray.)) +
    labs(x = "yearly dust flux [g/m-2/yr-1]", y = "difference in gypsum accumulation depth [meq/100g soil]",
    title = "sensitivity test for dust flux");

#sensitivity test for AESFactor
ETPlot = ggplot(data = results) + geom_smooth(se = FALSE, mapping = aes(x = results$AETFactor, y = results$difference.observedArray.)) +
    labs(x = "Evapotranspiration Factor", y = "difference in gypsum accumulation depth [meq/100g soil]",
    title = "sensitivity test for Evapotranspiration");

CalcGypsum(years = 10000, DustCa = 4, DustSO4 = 4, AETFactor = 200, Getresults = TRUE);
obsCalc


df = as.data.frame.array(obsCalc);
df$depth = (1:nrow(df) * 5);
sidebysidePlot =  ggplot(df) + geom_path(, mapping = aes(y = depth, x = Observed, colour = "Observed")) +
    geom_path(mapping = aes(y = depth, x = Calculated, colour = "Calculated")) + scale_y_reverse() + labs(x = "gypsum");


#find minimum value
minimal = temp[3, temp[1,] == min(unlist(temp[1,]))]

##create for Eilta
#EilatData = read.csv("Eilat.csv",);

##get amount values 
#EilatAmountCdf = ecdf(EilatData[2]$vals);

##Interval 
#EilatInterval = tail(EilatData$time, -1) - head(EilatData$time, -1);
#EilatInterval = EilatInterval[which(EilatInterval < 180)];

##get interval CDF
#EilatIntervalCDF = ecdf(EilatInterval);
#Intervalcdf = cbind(get("x", environment(EilatIntervalCDF)), get("y", environment(EilatIntervalCDF)));

##get CDF function
#Amountcdf = cbind(get("x", environment(EilatAmountCdf)), get("y", environment(EilatAmountCdf)));

##generate Rain series

#RainSeries = numeric(NumberOfYears * 365);
#interval = 0;
#RainSeries = sapply(RainSeries, FUN = function(X) GetRandomValue(X));

##create Module from CSM
#Rcpp::sourceCpp('C:/Users/Lior/master/evapocalc/Calcyp/main.cpp');
#n = new(CSM,2);


