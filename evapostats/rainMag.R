require(jsonlite)
require(httr)
require(lubridate)
stations = GET("https://api.ims.gov.il/v1/envista/stations", add_headers("Authorization" = "ApiToken f058958a-d8bd-47cc-95d7-7ecf98610e47"))
stations = (fromJSON(content(stations, as = "text")))
sdom = GET("https://api.ims.gov.il/v1/envista/stations/45/data/12/?from=2019/02/1&to=2020/05/03", add_headers("Authorization" = "ApiToken f058958a-d8bd-47cc-95d7-7ecf98610e47"))
sdom = (fromJSON(content(sdom, as = "text")))
tableSdomRaw = tibble(date = date(as.POSIXct(gsub("T", " ", sdom$data$datetime), tz = "GMT")), time = format(as.POSIXct(gsub("T", " ", sdom$data$datetime), tz = "GMT"), "%H:%M:%S"), value = sdom$data$channels %>% map_dbl(~.x$value)) %>% filter(value >= 0)
tableSdom = tableSdomRaw  %>% filter(value > 0.1, year(date) >= 2000) %>% rowwise() %>% mutate(intens = value * 6, runoffRatioS1 = GBRunoff(value * 6, 1), runoffRatioS2 = GBRunoff(value * 6, 2), runoffRatioS3 = GBRunoff(value * 6, 3)) %>% mutate(totRunoffS1 = runoffRatioS1 * value, totRunoffS2 = runoffRatioS2 * value, totRunoffS3 = runoffRatioS3 * value) %>% ungroup() %>% filter(intens<72)


elat = GET("https://api.ims.gov.il/v1/envista/stations/64/data/17?from=2002/11/29&to=2020/05/03", add_headers("Authorization" = "ApiToken f058958a-d8bd-47cc-95d7-7ecf98610e47"))
elat = (fromJSON(content(elat, as = "text")))
tableElatRaw = tibble(date = date(as.POSIXct(gsub("T", " ", elat$data$datetime))), time = format(as.POSIXct(gsub("T", " ", elat$data$datetime)), "%H:%M:%S"), value = elat$data$channels %>% map_dbl(~.x$value)) %>% filter(value >= 0)
tableElat = tableElatRaw %>% filter(value > 0.1) %>% rowwise() %>% mutate(intens = value * 6, runoffRatioS1 = GBRunoff(value * 6, 1), runoffRatioS2 = GBRunoff(value * 6, 2), runoffRatioS3 = GBRunoff(value * 6, 3)) %>% mutate(totRunoffS1 = runoffRatioS1 * value, totRunoffS2 = runoffRatioS2 * value, totRunoffS3 = runoffRatioS3 * value) %>% ungroup() #%>% filter(intens < 72)


daily = tableSdom %>% group_by(date) %>% summarise(daily = sum(value), Pleistocene.106k = sum(totRunoffS1), Pleistocene.31k = sum(totRunoffS2), Holocene.10k = sum(totRunoffS3), mean = mean(intens), max = max(intens), median = median(intens))

iqr = IQR(daily$daily)
Q <- quantile(daily$daily, probs = c(.25, .75), na.rm = FALSE)

daily = daily %>% filter(daily < Q[2] + 1.5 * iqr, daily > Q[1] - 1.5 * iqr)
lol = list(Holo = lm(Holocene.10k ~ 0 + daily, daily)[[1]], P31k = lm(Pleistocene.31k ~ 0 + daily, daily)[[1]], P106k = lm(Pleistocene.106k ~ 0 + daily, daily)[[1]])

%>% gather("surface", "value", - date, - median, - max, - daily, - mean) %>% ggplot(aes(daily, value, color = surface)) + geom_point()
#'%ni%' <- Negate('%in%')
#daily = daily %>% filter(value %ni% s$out)
formula = y~a*x^b
ggplot(daily, aes(x = daily, y = value, color = surface)) + geom_point() + stat_smooth(method = "nls", formula = 'y~a*x^b', method.args = list(start = c(a = 1, b = 1)), se = FALSE) + geom_text(x = 600, y = 1, label = power_eqn(daily), parse = TRUE) +
     labs(y = "daily runoff[mm]", x = "daily rain[mm]",title = "Elat 2002-2020")
ggplot(tableSdom %>% filter(value>0), aes(y=intens)) + geom_boxplot(outlier.color = "blue") 
save(tableElatRaw, tableSdomRaw, file = "rainIntensMeasured.RData")

simRain = tibble(age = seq(10000, 110000, by = 100), rain = sample(unique(SynthRainE %>% filter(rain>0) %>% pull(rain)), 1001, replace = T)) %>% mutate(runoff = ageRainfunc(age, rain)) 
ggplot(simRain, aes(x = age, y = rain, color = runoff, fill = runoff,z = runoff)) + geom_point() + scale_color_gradientn(colours = rainbow(5)) 

GBRunoff = function(intens, profile) {
   
    y=0
    if (profile == 1) {
        y = 0.00871 * intens
    } else if (profile == 2) {
        y = 0.0058 * intens
    } else {
        y =  0.0046 * intens
    }
    return(y)
}

power_eqn = function(df, start = list(a = 300, b = 1)) {
    m = nls(Discharge ~ a * Age ^ b, start = start, data = df);
    eq <- substitute(italic(y) == a ~ italic(x) ^ b,
               list(a = format(coef(m)[1], digits = 2),
                    b = format(coef(m)[2], digits = 2)))
    as.character(as.expression(eq));
}