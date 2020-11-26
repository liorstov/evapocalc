require(parallel)
load("synth.RData")
cl <- makeCluster(5, type = "PSOCK")

clusterExport(cl, c("T2.1", "T2.10", "T2.1", "opt.AETF", "opt.dust", "opt.FC", "opt.sulfate", "opt.WP", "Zel1", "Zel2"))

clusterEvalQ(cl, {
    require(tidyverse)
    require(Rcpp)
    require(ggplot2)
    require(reshape2)
    require(zoo)
    require(fitdistrplus)

    require(tictoc)
    Rcpp::sourceCpp('C:/Users/liorst/source/repos/evapocalc/Calcyp/CSM.cpp', verbose = TRUE, rebuild = 0);
    b <<- new(CSMCLASS);
    source("C:/Users/liorst/source/repos/evapocalc/evapostats/Functions.R", encoding = "Windows-1252")

    opt.AETF = 1.2;
    opt.WP = 0.013;
    opt.FC = 0.1;
    opt.sulfate = 10;
    opt.dust = 3;
    seq.AETF = seq(0.8, 1.2, by = 0.4) * opt.AETF %>% rep(60);
    seq.WP = seq(0.8, 1.2, by = 0.4) * opt.WP %>% rep(60);
    seq.FC = seq(0.8, 1.2, by = 0.4) * opt.FC %>% rep(60);
    seq.rainSeq = seq(0.8, 1.2, by = 0.4) * opt.sulfate %>% rep(60);
    seq.dustSeq = seq(0.8, 1.2, by = 0.4) * opt.dust %>% rep(60);

    RainArr = seq(10, 20, by = 1);
    DustArr = c(seq(1, 26, by = 5));
    repetition = 1:5
    rainDustArray = as.matrix(crossing(RainArr, DustArr, repetition))
})
results = list()
results$Terace12 = list();
results$Qa3= list();


days = seq(20,20);
annual = seq(30, 150, by = 20);
PET =( seq(1000, 2300, by = 100))
array = tibble(crossing(days, annual, PET), mean = annual / days) %>% filter(mean >= 3 & mean <= 11)


for (i in 1:length(PET)) {
     #eilat
    SynthRainEP = GenerateSeries(station = 347700, stationEvap = 347704,NumOfSeries = 1000, PETfactor = PET[i])


    #sedom
    SynthRainSP = GenerateSeries(station = 337000, stationEvap = 337000,NumOfSeries = 1000, PETfactor = PET[i])

  
    SynthRainEP %>% filter() %>% group_by(year) %>% summarise(s = sum(rain), n = n(), sum(PET)) %>% ungroup %>% summarise_all("mean")

    clusterExport(cl, "SynthRainEP")
    clusterExport(cl, "SynthRainSP")


    #rain sens ---

    results$Qa3 = c(results$Qa3, parLapply(cl, 1:30, fun = function(X) CalcGypsum(SynthRainEP,SynthRainEP, duration = 60000, Depth = 150)))
    results$Terace12 = c(results$Terace12, parLapply(cl, 1:30, fun = function(X) CalcGypsum(SynthRainSP,SynthRainSP, duration = 60000, Depth = 150)))
    save(results, file = "resultsListCalibRainEvapPleistocene.RData")

}


resultsTable = bind_rows(
   RMSD.T.10 = RectanglingResults(results$Qa3, T1.10) %>% mutate(profile = "T1.10", isHolocene = T),
    RMSD.Zel11 = RectanglingResults(results$Terace12, Zel11) %>% mutate(profile = "zel11", isHolocene = T),
#RMSD.T1.1 = RectanglingResults(results$T1.1, c(T1.1Observed %>% pull(correctedMean))) %>% mutate(profile = "T1.1", isHolocene = F),
#RMSD.T2.1 = RectanglingResults(results$T2.1, c(T2.1 %>% pull(correctedMean))) %>% mutate(profile = "T2.1", isHolocene = F),
#RMSD.zel1 = RectanglingResults(results$zel1, c(Zel1 %>% pull(correctedMean))) %>% mutate(profile = "zel1", isHolocene = F)
)
resultsTable = resultsTable %>% mutate(optimal = ifelse(FC == opt.FC & AETF == opt.AETF & sulfate == opt.sulfate & dustFlux == opt.dust & WP == opt.WP, T, F)) %>%
    mutate(region = ifelse(profile %in% c("zel11", "zel1"), "Zeelim", "Shehoret"))

## sens test for weather series---
#%>% group_by(scale = round(scale,1), shape = round(shape,1), region)%>% summarise_if(is.numeric,median) 
sensTest = resultsTable %>% filter(AnnualRain <=100) %>% mutate(meanDay = (AnnualRain / rainDays), Gypsum_depth = PeakDepth, Total_concentration = total) %>%
    gather("target", "value", Gypsum_depth, Total_concentration)

a = sensTest %>% filter(target != "Gypsum_depth") %>% ggplot(aes(x = AnnualRain, y = value, color = rainDays)) + geom_point(shape = 1, size = 5, show.legend = T) +facet_wrap(. ~ paste(region, ": ", duration, " years")) + labs() + coord_cartesian(xlim = c()) + labs(y = "gypsum conc.\n [meq/100 gr soil]", x = "Annual PET [mm]", color = "# rain days", size = "Rain event depth [mm]") + scale_color_distiller(palette = "PuBuGn", direction = 1) + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
b = sensTest %>% filter(target == "Gypsum_depth") %>% ggplot(aes(x = AnnualRain, y = value, color = rainDays)) + geom_point(shape = 1, size = 5) + facet_wrap(~region) + labs() + coord_cartesian(ylim = c()) + scale_y_reverse(name = "Depth [cm]") + xlab("Mean annual rainfall [mm]") + labs(size = "# rain days", color = "Rain event depth [mm]") + scale_color_distiller(palette = "PuBuGn", direction = 1) + theme(legend.position = "none",, strip.text = element_blank())
ggarrange(a,b,ncol =1)
sensTest %>% filter(target == "Gypsum_depth") %>% ggplot(aes(x = AnnualRain, y = value, color = factor(duration), group = (PET %/% 300 * 300))) + geom_boxplot() + geom_point(shape = 1, size = 5) + facet_wrap(~region) + labs() + coord_cartesian(ylim = c(150, 0)) + scale_y_reverse(name = "Depth [cm]") + xlab("Annual PET [mm]") + labs(color = "# rain days", size = "Rain event depth [mm]") # + scale_color_gradientn(colours = rainbow(5, rev = T))
a = sensTest %>% filter(target == "Gypsum_depth") %>% ggplot(aes(x = factor(PET %/% 300 * 300), y = value)) + geom_boxplot(outlier.color = NA) + stat_summary(geom = "point",fun = "mean", shape = 5) + facet_wrap(~region) + labs() + coord_cartesian(ylim = c()) + scale_y_reverse(name = "Depth [cm]") + xlab("Annual PET [mm]") + labs(color = "# rain days", size = "Rain event depth [mm]") # + scale_color_gradientn(colours = rainbow(5, rev = T))
sensTest %>% filter(target != "Gypsum_depth") %>% ggplot(aes(x = factor(PET %/% 300 * 300), y = value)) + geom_boxplot() + facet_wrap(~region) + labs() + coord_cartesian(xlim = c()) + labs(y = "mean gypsum conc. [meq/100 gr soil]", x = "Annual PET [mm]", color = "# rain days", size = "Rain event depth [mm]") # + scale_color_gradientn(colours = rainbow(5, rev = T))
b = sensTest %>% filter(target != "Gypsum_depth") %>% ggplot(aes(x = factor(PET %/% 300 * 300), y = value)) + geom_boxplot(outlier.color = NA) + stat_summary(geom = "point", fun = "mean", shape = 5) + facet_wrap(~region) + labs() + coord_cartesian(xlim = c()) + labs(y = "mean gypsum conc. [meq/100 gr soil]", x = "Annual PET [mm]", color = "# rain days", size = "Rain event depth [mm]") # + scale_color_gradientn(colours = rainbow(5, rev = T))
arrangeGrob(a,b)
grid.arrange(arrangeGrob(a, b))
#---
grid.arrange
