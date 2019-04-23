require(aqp)
require(RODBC)
require(lattice)
require(latticeExtra)
require(ggplot2)
require(ggpmisc)
setwd("C:\\Users\\Lior\\master\\evapocalc\\aqp")


without.RC.Horizons = TRUE
without.ParentMaterial.Sand = FALSE
without.numbers = TRUE
Ca.middle = FALSE

#declare functions
mean.and.sd <- function(values){
  m <- mean(values, na.rm=TRUE)
  s <- sd(values, na.rm=TRUE)
  upper <- m + s
  lower <- m - s
  res <- c(lower=lower, mean=m, upper=upper, sd=s)
  return(res)
}

#remove std dev from plots
remove.stdev = function(plot){
  plot$panel.args.common$upper = 0
  plot$panel.args.common$cf.col='white'
  return(plot)
}

#find the top most Bk horizon for each site
Get.top.ca.horizon <- function(horizons.table){
  #filter by ca horizons only
  dataframeCa = subset(horizons.table,subset =  grepl("(ca)|(k)",horizons.table$horizon, ignore.case=TRUE))
  dataframeCa = subset(dataframeCa,subset =  !grepl("[A]",dataframeCa$horizon, ignore.case=FALSE))
  dataframeCa = subset(dataframeCa,subset =  !grepl("[b]",dataframeCa$horizon, ignore.case=FALSE))
  
  #sort DF by roof
  dataframeCa = dataframeCa[order(dataframeCa$siteid,dataframeCa$horizon,dataframeCa$depthroof),]
  
  ca.horizons.count.per.site = aggregate(dataframeCa[,] , by=list(dataframeCa$siteid), length)[,1:2]
  colnames(ca.horizons.count.per.site)[2] = "count"
  top.horizon = aggregate(dataframeCa[,] , by=list(dataframeCa$siteid), head,1)
  
  #filter by horizons above 100 cm
  top.horizon = top.horizon[top.horizon$depthroof <= 100,]
  
  if (Ca.middle)
    #calculate middle of horizon
    top.horizon$middle_depth = (top.horizon$depthbase + top.horizon$depthroof)/2
  
  # aggregate by site
  return(merge(top.horizon,ca.horizons.count.per.site,by = "Group.1"))
}

get.sand.group = function(total.sand){
  if (is.na(total.sand)){
    return("Total sand <= 50%")
  }
  if (total.sand <= 50){
    return("Total sand <= 50%")
  }
  else{
    return("Total sand > 50%")
  }
}


##code 

plot.title = "Depth vs. median CaCO3 \n grouped by annual precipitation"

##connect to dns soil
remove(con)
con = odbcConnect("soils")

##create dataframe from horizons
dataframe = sqlQuery(con, "select horizons.siteid,horizons.horizon,horizons.depthroof,horizons.depthbase,horizons.caco3,horizons.caso4,horizons.Sieving_and_Pipette_Total_Sand,horizons.ColorMunsell,sites.MeanAnnualPrecipitation,sites.ParentTexture1,sites.SiteName,sites.ITM_X_coordinate,sites.ITM_Y_coordinate,sites.AgeClass1 from horizons,sites where sites.siteid=horizons.siteid and (caco3>=0 and caco3 is not null) and (horizons.siteid IN (select sites.siteid from sites where regionid = 94))")

if (without.numbers){
  dataframe = subset(dataframe,subset =  grepl("[a-z]",dataframe$horizon, ignore.case=TRUE))
  plot.title=paste(plot.title,"\nexclude numbers")
}

if (without.ParentMaterial.Sand){
  dataframe = subset(dataframe,(dataframe$ParentTexture1 != "sand" | is.na(dataframe$ParentTexture1)))
  plot.title=paste(plot.title,"\nexclude sand")
  }
if (without.RC.Horizons){
  dataframe = subset(dataframe,subset =  !grepl("[CRb]",dataframe$horizon, ignore.case=FALSE))
  plot.title=paste(plot.title,"\nexclude C/R/b")
}


#for some reason some rows are NA so i omit them 
#dataframe = na.omit(dataframe)

#divide to group according to percipation
range = c(0,80,120,200,250,300,350,400,450,500)
AP=as.numeric(dataframe$MeanAnnualPrecipitation)
dataframe$group = '1 - 80 mm'
dataframe$groupnumber = 1
dataframeLimit = length(dataframe$siteid)
for (i in 1:dataframeLimit){
  for (j in 2:length(range)){
    #validate percipitation availability
    if (!is.na(dataframe$MeanAnnualPrecipitation[i]) & dataframe$MeanAnnualPrecipitation[i] <= range[j]){
      dataframe$group[i] = paste(range[j-1]+1, '-',range[j], 'mm')
      dataframe$groupnumber[i]=j-1
      j=length(range)
      break
    }
    
  }
}

dataframe$group[which(dataframe$group == '1 - 80 mm')] = '< 80 mm'



#set 0 sand to NA
dataframe$Sieving_Pipette_Total_Sand[which(dataframe$Sieving_Pipette_Total_Sand==0)]=NA

#add column for sand groups
dataframe$sand.group = mapply(FUN = get.sand.group, dataframe$SievingPipetteTotalSand)

sumColums=table(dataframe$siteid,dataframe$groupnumber)
sumColums =colSums(sumColums != 0)



group.list=unique(dataframe$group[order(dataframe$groupnumber)])

#get the shallowest ca horizon for each site
top.ca.per.site = Get.top.ca.horizon(dataframe)
group.ca.count = colSums(table(top.ca.per.site$siteid, top.ca.per.site$groupnumber) != 0)


#get a list of all the sites
site.list = aggregate(dataframe[,] , by=list(dataframe$siteid), head,1)

#promote to SoilProfileCollection
temp = dataframe
depths(dataframe)<- siteid ~ depthroof + depthbase
site(dataframe) <- ~groupnumber

#slub.structure = horizon thickness cm

caso4.slab <- slab(dataframe, fm = groupnumber ~ caso4, slab.structure = 1, strict = FALSE)
caco3.slab <- slab(dataframe, fm = groupnumber ~ caco3, slab.structure = 1, strict = FALSE)
cac03.slab.siteid = slab(dataframe, fm = siteid ~ caco3, slab.structure = 1, strict = FALSE, slab.fun = mean.and.sd)
sand.slab <- slab(dataframe, fm = groupnumber ~ SievingPipetteTotalSand, slab.structure = 1, strict = FALSE, slab.fun = mean.and.sd)
t = merge(caso4.slab, temp[, c("siteid", "groupnumber")], by = "siteid")

test=mapply(FUN = function(x,y,z) paste0(x,'\n','CaCO3 Sites = ', y, '\nBk sites = ',z),group.list,as.character(sumColums), as.character(group.ca.count))

#create plot
my.plot <- xyplot(top ~ p.q50 | as.character(groupnumber), data = caso4.slab, main = list(label = plot.title, cex = 0.75), ylab = 'Depth [cm]', xlab = 'Median CaCO3 concentration bounded by standard deviation [%]',
                  lower = caso4.slab$p.q25, upper = caso4.slab$p.q75, sync.colors = TRUE, alpha = 0.25,
                  cf = caso4.slab$contributing_fraction, cf.col = 'black',
                  ylim=c(105,-5), xlim=c(0,70), layout=c(3,1),
                  par.settings=list(layout.heights = list(strip=2.30),superpose.line=list(lwd=3.5, col=c('blue')),layout.widths = list(6)),
                  panel=panel.depth_function,
                  prepanel=prepanel.depth_function,
                  par.strip.text = list(cex=0.8, lines = 1.6),
                  auto.key=list(columns=5, lines=TRUE, points=FALSE), strip = strip.custom(factor.levels = paste(test)))

my.plot1 <- xyplot(top ~ mean , group = groupnumber, data=caco3.slab, main = list(label = plot.title, cex=0.75),ylab = 'Depth [cm]', xlab = 'Mean CaCO3 concentration bounded by standard deviation [%]',
                  lower=caco3.slab$lower, upper=caco3.slab$upper, sync.colors=TRUE, alpha=0.25,
                  cf=caco3.slab$contributing_fraction, cf.col = 'black',
                  ylim=c(105,-5), xlim=c(0,70), layout=c(1,1),type=c('S','g'),
                  par.settings=list(layout.heights = list(strip=2.30),superpose.line=list(lwd=2),layout.widths = list(6)),
                  auto.key=list(columns=5, lines=TRUE, points=FALSE), strip = strip.custom(factor.levels = paste(test)))

sand.plot <- xyplot(top ~ mean | as.character(groupnumber), data=sand.slab, main = list(label = "Depth vs. mean total sand \n grouped by winter precipitation percentage \n", cex=0.75),ylab = 'Depth [cm]', xlab = 'Mean sand percentage bounded by standard deviation [%]',
                  lower=sand.slab$lower, upper=sand.slab$upper, sync.colors=TRUE, alpha=0.25,
                  cf=sand.slab$contributing_fraction, cf.col = 'black',
                  ylim=c(105,-5), xlim=c(0,100), layout=c(3,1), scales = list(x=list(tick.number = 4)),
                  par.settings=list(layout.heights = list(strip=2.30),superpose.line=list(lwd=2, col=c('red', 'Orange2')),layout.widths = list(6)),
                  panel=panel.depth_function,
                  prepanel=prepanel.depth_function,
                  par.strip.text = list(cex=0.8, lines = 1.6),
                  auto.key=list(columns=5, lines=TRUE, points=FALSE), strip = strip.custom(factor.levels = paste(test)))


top.ca.plot <- xyplot((depthroof ~ caco3 | as.character(groupnumber)), data = top.ca.per.site, col = c("red"),
                      panel=function(x,y,...) { 
                        panel.xyplot(x,y,...) 
                        panel.abline(h=median(y))} )


temp.plot = xyplot(top ~ mean | as.character(groupnumber), data=t, groups = siteid, type =c("l") ,ylim=c(105,-5), xlim=c(0,70),
                   auto.key=list(columns=5, lines=TRUE, points=FALSE), strip = strip.custom(factor.levels = paste(test)),
                   par.settings = list(superpose.line = list(col=c('royal blue'), lwd = 2)),
                   panel=function(x, y, ...) {
                     panel.xyplot(x, y, ...);
                     ltext(x=5, y=10, labels=as.integer(NULL), pos=1, offset=1, cex=0.8  )})

formula <- y ~ poly(x, 2, raw = TRUE)
top.ca.to.rain <- ggplot(data = top.ca.per.site, aes(y=depthroof, x=MeanAnnualPrecipitation, color = sand.group),)
top.ca.to.rain = top.ca.to.rain + ggtitle(label = paste("Top Ca horizons vs. mean Annual Precipitation \n",c(length(top.ca.per.site[,1])),"Ca sites")) + xlab("Mean annual precipitation [mm/y]")+ylab("Depth [cm]")+geom_point() + scale_y_reverse() +  geom_smooth(se=FALSE, method = "lm", formula = formula) +
  stat_poly_eq(
    aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
    formula = formula, parse = TRUE,  label.y.npc  = 0.1
  ) 
  
my.plot = my.plot + as.layer(top.ca.plot) +as.layer(temp.plot)

double.plot = remove.stdev(my.plot)+as.layer(remove.stdev(sand.plot))


write.csv(x = top.ca.per.site,file = "SlabOutput.csv")
pdf(file = paste(format(Sys.time(), "%b_%d_%Y"),"soilPlot.pdf"))
print(my.plot)
print(sand.plot)
print(top.ca.to.rain)
print(double.plot)
dev.off()

my.plot
