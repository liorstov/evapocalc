require(lubridate)
require(tibble)
require(MASS)
require(actuar)

setwd('C:/Users/Lior/master/evapocalc/evapostats');
IMSRain = as_tibble(read.csv("Eilat.csv")[seq(1, 4)]);

pd = fitdistr(IMSRain$vals, "weibull")

IMSRain$serial = apply(IMSRain,1, function(X) {
    lubridate::yday(as.POSIXct.Date(X[1] - 719529, origin = '1970-01-01'));
})
IMSRain$year = apply(IMSRain, 1, function(X) {
    lubridate::year(as.POSIXct.Date(X[1] - 719529, origin = '1970-01-01'));
})

numOfYears = length(unique(IMSRain$year));

IMSRain$IsWAW = ((dplyr::lag(IMSRain$time) - IMSRain$time) == -1)

pWAW = hist(IMSRain$serial[which(IMSRain$IsWAW)], 0:365)
pWAD = hist(IMSRain$serial[which(!IMSRain$IsWAW)], 0:365)
pWet = hist(IMSRain$serial, 0:365)

ProbRainDays = tibble(DayNum = seq(1:365), WAWCount = pWAW$counts,
                        WADCount = pWAD$counts,
                        WetCount = pWet$counts);


ProbRainDays$PrevWet = dplyr::lag(ProbRainDays$WetCount)

ProbRainDays$pWAW = ProbRainDays$WAWCount / ProbRainDays$PrevWet;
ProbRainDays$pWAD = ProbRainDays$WADCount / (numOfYears - ProbRainDays$PrevWet);

ProbRainDays$pWAW = 0;
ProbRainDays$pWAD = 0;

dinvweibull(1, shape = pd$estimate[1], scale = pd$estimate[2])

DepthLimit = 0.1;

SyntRain = tibble(Prob = runif(365000), Depth = NA)

if (SyntRain$Prob[1] <= ProbRainDays$pWAD[1]) {
    SyntRain$Depth = dweibull(runif(1), pd$estimate[1], pd$estimate[2]) + DepthLimit;
} else {
    SyntRain$Depth = 0;
}

for (I in seq(1:) {

}