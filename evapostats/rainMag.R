require(jsonlite)
require(httr)
require(lubridate)
stations = GET("https://api.ims.gov.il/v1/envista/stations", add_headers("Authorization" = "ApiToken f058958a-d8bd-47cc-95d7-7ecf98610e47"))
stations =(fromJSON(content(msg1, as = "text")))
sdom = GET("https://api.ims.gov.il/v1/envista/stations/65/data/12/?from=2001/11/29&to=2020/05/03", add_headers("Authorization" = "ApiToken f058958a-d8bd-47cc-95d7-7ecf98610e47"))
sdom = (fromJSON(content(sdom, as = "text")))
tableSdom = tibble(date = date(as.POSIXct(gsub("T", " ", sdom$data$datetime))), time = format(as.POSIXct(gsub("T", " ", sdom$data$datetime)),"%H:%M:%S"), value = sdom$data$channels %>% map_dbl(~.x$value))



elat = GET("https://api.ims.gov.il/v1/envista/stations/64/data/17?from=2002/11/29&to=2020/11/30", add_headers("Authorization" = "ApiToken f058958a-d8bd-47cc-95d7-7ecf98610e47"))
elat = (fromJSON(content(elat, as = "text")))
tableElat = tibble(date = date(as.POSIXct(gsub("T", " ", elat$data$datetime))), time = format(as.POSIXct(gsub("T", " ", elat$data$datetime)), "%H:%M:%S"), value = elat$data$channels %>% map_dbl(~.x$value))
table = table %>% filter(value >= 0)
daily = table %>% mutate(magn = value * 6) %>% group_by(date) %>% summarise(daily = sum(value), mean = mean(magn), max = max(magn), median = median(magn))%>%dplyr::select(-date) %>% gather("param","value",-daily)
ggplot(daily %>% filter(daily > 0), aes(x = daily, y = value, color = param)) + geom_point() + coord_cartesian(ylim = c(0, 3))

save("mag.RData")