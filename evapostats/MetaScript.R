
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


#set wd
setwd("C:/Users/Lior/master/evapocalc/evapostats")

#init variables
NumberOfYears = 100;


#create for Eilta
EilatData = read.csv("Eilat.csv",);

#get amount values 
EilatAmountCdf = ecdf(EilatData[2]$vals);

#Interval 
EilatInterval = tail(EilatData$time, -1) - head(EilatData$time, -1);
EilatInterval = EilatInterval[which(EilatInterval < 180)];

#get interval CDF
EilatIntervalCDF = ecdf(EilatInterval);
Intervalcdf = cbind(get("x", environment(EilatIntervalCDF)), get("y", environment(EilatIntervalCDF)));

#get CDF function
Amountcdf = cbind(get("x", environment(EilatAmountCdf)), get("y", environment(EilatAmountCdf)));

#generate Rain series

RainSeries = numeric(NumberOfYears * 365);
interval = 0;
RainSeries = sapply(RainSeries, FUN = function(X) GetRandomValue(X));


