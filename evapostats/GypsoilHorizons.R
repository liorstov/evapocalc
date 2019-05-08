require(aqp)
require(RODBC)
require(lattice)
require(latticeExtra)
require(ggplot2)
require(ggpmisc)
require(dplyr)
setwd("C:\\Users\\Lior\\master\\evapocalc\\aqp")


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


##code 

plot.title = "CaSO4 concentration for various profiles"

##connect to dns soil
remove(con)
con = odbcConnect("soils")

##create dataframe from horizons
dataframe = tbl_df(sqlQuery(con, "select horizons.siteid,horizons.horizon,horizons.depthroof,horizons.depthbase,horizons.caco3,horizons.caso4,horizons.Sieving_and_Pipette_Total_Sand,horizons.ColorMunsell,sites.MeanAnnualPrecipitation,sites.SiteName,sites.ITM_X_coordinate,sites.ITM_Y_coordinate,sites.AgeClass1,sites.AvgAge  from horizons,sites where (sites.siteid=horizons.siteid) and  (horizons.siteid IN (select sites.siteid from sites where regionid = 94))"));
ages = tbl_df(sqlQuery(con, "select ageclass,min_age,max_age from LUT_Age"));
dataframe = dataframe %>% full_join(ages, by = c("AgeClass1" = "ageclass"))
#omit na
dataframe = dataframe[which(!is.na(dataframe$caso4)),];
dataframe = dataframe[which(!is.na(dataframe$AgeClass1)),];


# filter by arid environment
dataframe = dataframe[which(dataframe$MeanAnnualPrecipitation <= 100),]
#for some reason some rows are NA so i omit them 
#dataframe = na.omit(dataframe)

#divide to group according to percipation
range = unique(dataframe$MeanAnnualPrecipitation);
dataframe$group = apply(dataframe[, c("MeanAnnualPrecipitation")], 1, FUN = function(X) getPercipitationRange(X, range));
dataframe$groupNum = as.integer(factor(dataframe$group))


##get the shallowest ca horizon for each site
#top.ca.per.site = get.top.ca.horizon(dataframe)
#group.ca.count = colsums(table(top.ca.per.site$siteid, top.ca.per.site$groupnumber) != 0)


#get a list of all the sites
site.list = aggregate(dataframe[,] , by=list(dataframe$siteid), head,1)

#promote to SoilProfileCollection

AQPClass = as.data.frame(dataframe);
depths(AQPClass) <- siteid ~ depthroof + depthbase 
site(AQPClass) <- ~ AgeClass1

#slub.structure = horizon thickness cm

caso4.slab <- slab(AQPClass, fm = groupNum ~ caso4, slab.structure = 1, strict = FALSE)
caso4.slab.siteid = slab(AQPClass, fm = AgeClass1 ~ caso4, slab.structure = 1, strict = FALSE, slab.fun = mean.and.sd)
caso4.slab.siteid = merge(caso4.slab.siteid, temp[, c("siteid", "MeanAnnualPrecipitation")], by = "siteid")
#this is for the strips
s
#create plot

my.plot1 = xyplot(top ~ mean | paste("siteid: ",siteid, "\nMAP: ", MeanAnnualPrecipitation), data = caso4.slab.siteid, lower = caso4.slab.siteid$lower, upper = caso4.slab.siteid$upper, main = list(label = plot.title, cex = 0.75), ylab = 'Depth [cm]', xlab = 'CaSO4 concentration [meq/100g soil]',
                        ylim = c(105, -5), xlim = c(-10, 70), layout = c(4, 1),
                        panel = panel.depth_function,
                  prepanel = prepanel.depth_function,
            par.strip.text = list(cex = 0.6, lines = 2.2),
                   auto.key = list(columns = 5, lines = TRUE, points = FALSE), strip = strip.custom())

my.plot2 = xyplot(top ~ mean | paste("MAP: ", MeanAnnualPrecipitation), data = caso4.slab.siteid, lower = caso4.slab.siteid$lower, upper = caso4.slab.siteid$upper, main = list(label = plot.title, cex = 0.75), ylab = 'Depth [cm]', xlab = 'CaSO4 concentration [meq/100g soil]',
                        ylim = c(105, -5), xlim = c(-10, 70), layout = c(4, 1),
                        panel = panel.depth_function,
                  prepanel = prepanel.depth_function,            
                   auto.key = list(columns = 5, lines = TRUE, points = FALSE), strip = strip.custom())

my.plot3 = xyplot(top ~ mean | paste("Age: ", AgeClass1), data = caso4.slab.siteid, lower = caso4.slab.siteid$lower, upper = caso4.slab.siteid$upper, main = list(label = plot.title, cex = 0.75), ylab = 'Depth [cm]', xlab = 'CaSO4 concentration [meq/100g soil]',
                        ylim = c(105, -5), xlim = c(-10, 70), layout = c(4, 1),
                        panel = panel.depth_function,
                  prepanel = prepanel.depth_function,
            par.strip.text = list(cex = 0.6, lines = 2.2),

                   auto.key = list(columns = 5, lines = TRUE, points = FALSE), strip = strip.custom())

ggplot(data = dataframe, mapping = aes(x = AgeClass1, y = MeanAnnualPrecipitation, group = siteid)) + geom_point()

#write.csv(x = top.ca.per.site,file = "SlabOutput.csv")
pdf(file = paste(format(Sys.time(), "%b_%d_%Y"),"soilPlotage.pdf"))
print(my.plot3)
dev.off()

