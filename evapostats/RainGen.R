require(lubridate)
require(tibble)
require(MASS)
require(actuar)
library(dplyr)



remove(con)
con = odbcConnect("RainEvap")

setwd('C:/Users/Lior/master/evapocalc/evapostats');

#get percipitation for 1988 to current
IMSRain = tbl_df(sqlQuery(con, "SELECT *, year(time) as year,month(time) as month ,dayofyear(time) as dayIndex  FROM data_dream.precip_daily where (precip_daily.idstation = 347700 and year(precip_daily.time) >= 1988)"));

#only rain days
IMSRain = IMSRain[IMSRain$measure > 0,]

pd = fitdistr(IMSRain$measure, "weibull")

#get wet after wet
IMSRain$prevDay = IMSRain$dayIndex - lag(IMSRain$dayIndex, default = 0);
IMSRain$WAW = IMSRain$prevDay == 1 | IMSRain$prevDay == -364

#frequency of wet days
WAWProb =  IMSRain %>% filter(WAW) %>% group_by(dayIndex) %>% summarise(WAWP = n()/nrow(IMSRain))
WADProb = IMSRain %>% filter(!WAW) %>% group_by(dayIndex) %>% summarise(WADP = n() / nrow(IMSRain))
WtDaysProb = WAWProb %>% full_join(WADProb)
WtDaysProb = WtDaysProb %>% replace_na(list(WAWP = 0, WADP = 0))

fitdistr(WtDaysProb$WAWP, "normal")

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

