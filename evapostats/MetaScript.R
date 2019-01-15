
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

CalcGypsum <- function(raindata = EilatData[, 1], years = 100, Depth = 50, thick = 5, wieltingPoint = 0.039, InitialCa = 10, initialSO4 = 10
                       , BulkDensity = 1.44, nArea = 1, FieldCapacity = 0.2, DustCa = 0, DustSO4 = 0, observedArray = observedCom) {
    b = new(CSMCLASS);
    list = b$Calculate(raindata, years, Depth, thick, wieltingPoint, InitialCa, initialSO4, BulkDensity, nArea, FieldCapacity,
                   DustCa, DustSO4);
    #cat(wieltingPoint, "      ", list$gypsum, '\n')
    observedArray[, 2] = list$gypsum;
    rm(list);
    return(list(RMSD = rmsd(observedArray), Result = observedArray));
}

require(Rcpp)
require(ggplot2)
require(reshape2)

#set wd
setwd("C:/Users/Lior/master/evapocalc/");
EilatData = read.csv("DB/rainseriestemp.csv");
Observed = read.csv("DB/measured.CSV");






#create observed comartments
observedCom = matrix(nrow=10, ncol = 2);
observedCom[1,1] = 0.72;
observedCom[2:8,1] = 2.44;
observedCom[8:10,1] = 4.2;

Rcpp::sourceCpp('C:/Users/Lior/master/evapocalc/Calcyp/CSM.cpp', verbose = TRUE);
temp = CalcGypsum(years = 100);
wiltingPointArray = seq(0.001, 0.04, 0.001);


temp = lapply(wiltingPointArray, FUN = function(X) CalcGypsum( wieltingPoint = X));

#find minimum value
minimal = temp[2, temp[1,] == min(unlist(temp[1,]))]

plot(x= wiltingPointArray, y = temp);
df = as.data.frame.array(minimal$Result);
df$depth = (1:nrow(df) * thick);
ggplot(df) + geom_path(, mapping = aes(y = depth, x = Observed, colour = "Observed")) +
    geom_path(mapping = aes(y = depth, x = Calculated, colour = "Calculated")) + scale_y_reverse() + labs(x = "gypsum");



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


