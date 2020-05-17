require(aqp)
require(RODBC)
require(lattice)
require(latticeExtra)
require(ggplot2)
require(ggpmisc)
require(dplyr)
require(tidyverse)

setwd("C:\\Users\\liorst\\Source\\Repos\\evapocalc")


without.RC.Horizons = FALSE
without.ParentMaterial.Sand = FALSE
without.numbers = FALSE
Ca.middle = FALSE

#declare functions
mean.and.sd <- function(values) {    
  m <- mean(values, na.rm=TRUE)
  s <- sd(values, na.rm=TRUE)
  upper <- m + s
  lower <- m - s
  res <- c(lower=lower, mean=m, upper=upper, sd=s)
  return(res)
}


getPercipitationRange <- function(MAP, range) {
    range = sort(range);
    Upper = head(which(range >= MAP), n = 1);
    Lower = Upper - 1;    

    return(paste(range[Lower], '-', range[Upper], 'mm'));
}

setAge <- function(class, avg, min, max) { 
    if (is.na(avg) & !is.na(class)) {
        avg = (as.numeric(min) + as.numeric(max)) %/% 2;
    }
    return(as.numeric(avg));
}


##code 

plot.title = "CaSO4 concentration for various profiles"

##connect to dns soil
remove(con)
con = odbcConnect("soils")

##create dataframe from horizons
dataframe = as_tibble(sqlQuery(con, "select horizons.siteid,horizons.horizon,horizons.depthroof,horizons.depthbase,horizons.caco3,horizons.caso4,horizons.Sieving_and_Pipette_Total_Sand,horizons.ColorMunsell,sites.MeanAnnualPrecipitation,sites.SiteName,sites.ITM_X_coordinate,sites.ITM_Y_coordinate,sites.AgeClass1 as ageclass,sites.AvgAge,sites.SoilType, list_desert.DesertID,list_desert.subregion  from horizons,sites,list_desert where (sites.siteid=horizons.siteid) and (sites.DesertID=list_desert.DesertID) and (horizons.siteid IN (select sites.siteid from sites where regionid = 94))"));
ages = as_tibble(sqlQuery(con, "select ageclass,min_age,max_age from LUT_Age"));
dataframe = dataframe %>% full_join(ages,"ageclass")
dataframe$AvgAge = apply(dataframe[c("ageclass", "AvgAge","min_age","max_age")], 1, function(X) setAge(X[1],X[2],X[3],X[4]))

#filter by regsols
dataframe = dataframe[which(dataframe$SoilType == 859),]

#divide to group according to percipation
range = unique(dataframe$MeanAnnualPrecipitation);
dataframe$group = apply(dataframe[, c("MeanAnnualPrecipitation")], 1, FUN = function(X) getPercipitationRange(X, range));
dataframe$groupNum = as.integer(factor(dataframe$group))

dataframe = dataframe %>% filter(!str_detect(SiteName, "12|13|14|15")) 


#get a list of all the sites
site.list = aggregate(dataframe[,] , by=list(dataframe$siteid), head,1)

#promote to SoilProfileCollection

AQPClass = as.data.frame(dataframe);
depths(AQPClass) <- SiteName ~ depthroof + depthbase 
site(AQPClass) <- ~ subregion

#slub.structure = horizon thickness cm

#caso4.slab <- slab(AQPClass, fm = groupNum ~ caso4, slab.structure = 1, strict = FALSE)
caso4.slab.siteid = as_tibble(slab(AQPClass, fm = SiteName ~ caso4, slab.structure = 5, strict = FALSE, slab.fun = mean.and.sd))
#caso4.slab.siteid = slab(AQPClass, fm = subregion ~ caso4, slab.structure = 1, strict = FALSE, slab.fun = mean.and.sd)


caso4.slab.siteid = caso4.slab.siteid  %>% left_join(dataframe[, c("SiteName", "AvgAge", "ageclass", "MeanAnnualPrecipitation")], by = "SiteName")
#this is for the strips

#data for plotting 
dataplot = caso4.slab.siteid[grep("ZEL", caso4.slab.siteid$SiteName),]
dataplot = dataplot[grep("early", dataplot$ageclass),]
dataplot = (caso4.slab.siteid)
#create plot

my.plot1 = xyplot(top ~ mean | paste(SiteName, "\n", ageclass), data = dataplot, upper = caso4.slab.siteid$upper, ylab = list('Depth [cm]', cex = 3), xlab = list('CaSO4 concentration [meq/100g soil]', cex = 3),
                        ylim = c(105, -5), layout = c(),
                        panel = panel.depth_function, cex = 9,
                  prepanel = prepanel.depth_function,
            par.strip.text = list(cex = 2.2, lines = 3),
            scales = list(x = list(tick.number = 2, cex = 2), y = list(cex = 2)),
            par.settings = list(superpose.line = list(col = 'RoyalBlue', lwd =5)),
                   auto.key = list(lines = TRUE, points = FALSE), strip = strip.custom(),
                   index.cond = list(c(1,3,7,4,6,2,5,8)))

my.plot2 = xyplot(top ~ mean | paste("MAP: ", MeanAnnualPrecipitation), data = caso4.slab.siteid, lower = caso4.slab.siteid$lower, upper = caso4.slab.siteid$upper, main = list(label = plot.title, cex = 0.75), ylab = 'Depth [cm]', xlab = 'CaSO4 concentration [meq/100g soil]',
                        ylim = c(105, -5), xlim = c(-10, 70), layout = c(4, 1),
                        panel = panel.depth_function,
                  prepanel = prepanel.depth_function,            
                   auto.key = list(columns = 5, lines = TRUE, points = FALSE), strip = strip.custom())

my.plot3 = xyplot(top ~ mean | paste(subregion), data = caso4.slab.siteid, lower = caso4.slab.siteid$lower, upper = caso4.slab.siteid$upper, main = "Gypsum aggregation", ylab = 'Depth [cm]', xlab = 'CaSO4 concentration [meq/100g soil]',
                        ylim = c(105, -5), layout=c(1,3),
                        panel = panel.depth_function,
                  prepanel = prepanel.depth_function,
            par.strip.text = list(cex = 0.6, lines = 2.2),
            index.cond = list(c(1,2)))

                   auto.key = list(columns = 5, lines = TRUE, points = FALSE), strip = strip.custom())

ggplot(data = dataframe, mapping = aes(x = (AvgAge), y = MeanAnnualPrecipitation, group = siteid)) + geom_point() +
  scale_x_continuous(breaks = seq(0,100000,5000) , limits = c(1,80000))

#write.csv(x = top.ca.per.site,file = "SlabOutput.csv")
pdf(file = paste(format(Sys.time(), "%b_%d_%Y"),"soilPlotage.pdf"),paper = "a4" )
print(my.plot3)  
dev.off()

