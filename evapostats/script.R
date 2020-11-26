
require(jsonlite)
require(httr)
require(lubridate)
require(tidyr)
require(tibble)


ziqron = GET("https://api.ims.gov.il/v1/envista/stations/45/data/daily/2015/02/04", add_headers("Authorization" = "ApiToken TOKEN_HERE"))
ziqron = (fromJSON(content(ziqron, as = "text")))

ziqron = tibble(date = rep(date(as.POSIXct(gsub("T", " ", ziqron$data$datetime), tz = "GMT"))), time = (format(as.POSIXct(gsub("T", " ", ziqron$data$datetime), tz = "GMT"), "%H:%M:%S")), value = ziqron$data$channels)
bla = ziqron %>% arrange(date, time) %>% unnest(value) %>% pivot_wider(c(date, time), names_from = name, values_from = value)
 
#combine different days
test = bind_rows(test,bla)
