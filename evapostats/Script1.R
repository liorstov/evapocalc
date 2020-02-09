plotLimit = 365;
histBreaksSize = 1;
#hist for observed
IMSRainDay = elat[[2]] %>% mutate(wet = (rain > 0) * 1) %>% group_by(dayIndex) %>% summarise(elat = sum(wet) / measuredYears)
#hist for simulated
IMSRainDayS = sedom[[2]] %>% mutate(wet = (rain > 0) * 1) %>% group_by(dayIndex) %>% summarise(sedom = sum(wet) / measuredYears)
SimRainDay = elat[[1]] %>% mutate(wet = (rain > 0) * 1) %>% group_by(dayIndex, SeriesNumber) %>% summarise(density = sum(wet) / measuredYears)
SimRainDay = SimRainDay %>% group_by(dayIndex) %>%
                summarise(Simulated = mean(density), min = quantile(density, 0.05), median = quantile(density, 0.5), max = quantile(density, 0.95)) %>%
                    add_column(month = (dmy("1-09-2000") + 0:364)) %>% left_join(IMSRainDay, by = "dayIndex") %>%
                       left_join(IMSRainDayS, by = "dayIndex");

SimRainDayS = sdom[[1]] %>% mutate(wet = (rain > 0) * 1) %>% group_by(dayIndex, SeriesNumber) %>% summarise(density = sum(wet) / measuredYears)
SimRainDayS = SimRainDayS %>% group_by(dayIndex) %>%
                summarise(Simulated = mean(density), min = quantile(density, 0.05), median = quantile(density, 0.5), max = quantile(density, 0.95)) %>%
                    add_column(month = (dmy("1-09-2000") + 0:364))

p1 = ggplot(data = SimRainDay, aes(x = month, group = 1)) + geom_line(aes(y = elat, color = "Elat")) + ylab("Wet Probability [-]") + xlab("") +
                                geom_line(aes(y = sedom, color = "Sedom"), size = 1.5) +
                                   scale_x_date(date_labels = "%B", expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0.0));

rainProbabilities = ggplot(data = SimRainDay, aes(x = month, group = 1)) +
                                geom_line(aes(y = PWET, color = "WET")) +
                                geom_line(aes(y = PWAW, color = "WAW")) +
                                geom_line(aes(y = PWAD, color = "WAD")) +
                                   scale_x_date(date_labels = "%B", expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0.0)) + ylab("Wet Probability [-]") + xlab("");



plotLimit = 365;
histBreaksSize = 1;
#hist for observed
IMSPetE = elat[[2]] %>% filter(!is.na(pen)) %>% group_by(dayIndex) %>% summarise(elat = mean(pen))
petYears = elat[[2]] %>% filter(!is.na(pen)) %>% distinct(waterYear) %>% nrow()
elat[[1]]$PETSeries = elat[[1]]$year %/% petYears;


IMSPetS = sedom[[2]] %>% filter(!is.na(pen)) %>% group_by(dayIndex) %>% summarise(sdom = mean(pen))
petYears = sedom[[2]] %>% filter(!is.na(pen)) %>% distinct(waterYear) %>% nrow()
sedom[[1]]$PETSeries = sedom[[1]]$year %/% petYears;



#aggregating   
SimPETDay = IMSPetE %>% left_join(IMSPetS, by = "dayIndex") %>% mutate(month = (dmy("1-09-2000") + dayIndex))
SimPETDay = sedom[[1]] %>% group_by(dayIndex) %>% summarise(Simulated = mean(PET, na.rm = TRUE), min = quantile(PET, 0.05, na.rm = TRUE), median = quantile(PET, 0.5, na.rm = TRUE), max = quantile(PET, 0.95, na.rm = TRUE)) %>%
                left_join(IMSPet, by = "dayIndex") %>% left_join(sedom[[4]], by = "dayIndex") %>% mutate(month = (dmy("1-09-2000") + dayIndex)) %>%
                mutate(K = replace(K, K == 1, "Wet"), K = replace(K, K == 2, "DAW"), K = replace(K, K == 3, "DAD"))
p4 = ggplot(data = SimPETDay, aes(x = month, group = 1)) +
            geom_line(aes(y = elat, color = "Eilat")) + ylab("PET [mm / day]") + xlab("") +
            geom_line(aes(y = sedom, color = "Sedom"), size = 1.5) +
             scale_x_date(date_labels = "%B", expand = c(0, 0), label = "") + scale_y_continuous(expand = c(0, 0.0)) 



plotLimit = 100;
histBreaksSize = 5;

#histogram for observed
IMSRainAnnE = elat[[2]] %>% filter(rain > 0.1) %>% group_by(waterYear) %>% summarise(annual = sum(rain), WetDays = n()) %>% filter(annual <= plotLimit)
IMShistannualE = hist(IMSRainAnnE$annual, breaks = seq(0, plotLimit, histBreaksSize), plot = 0) %>%
                        Hist2tibble() %>% mutate(obsDensityE = counts / measuredYears)


#histogram for observed
IMSRainAnnS = sedom[[2]] %>% filter(rain > 0.1) %>% group_by(waterYear) %>% summarise(annual = sum(rain), WetDays = n()) %>% filter(annual <= plotLimit)
IMShistannualS = hist(IMSRainAnnS$annual, breaks = seq(0, plotLimit, histBreaksSize), plot = 0) %>%
                        Hist2tibble() %>% mutate(obsDensityS = counts / measuredYears)
#histogram for simulated
ggplot(IMShistannualE %>% left_join(IMShistannualS, by = "breaks"), aes(x = breaks)) +
        ylab("pdf") +
        geom_line(aes(y = obsDensityE, color = "Elat")) +
        geom_line(aes(y = obsDensityS, color = "Sedom")) +
        xlab("Annual rain mm/year") +
        scale_y_continuous(expand = c(0, 0.0)) + scale_x_continuous(expand = c(0, 0));



#annual wet days----
plotLimit = 40;
histBreaksSize = 1;
#histogram for observed
IMSAnnWetDaysE = hist(IMSRainAnnE$WetDays, breaks = seq(-2, plotLimit, histBreaksSize), plot = 0) %>%
                        Hist2tibble() %>% rename(elat = density)
IMSAnnWetDaysS = hist(IMSRainAnnS$WetDays, breaks = seq(-2, plotLimit, histBreaksSize), plot = 0) %>%
                        Hist2tibble() %>% rename(sdom = density)


#histogram for simulated
#SimAnnWetDays = simRainAnn %>% group_by(SeriesNumber) %>%
                       #summarise(a = list(hist(WetDays, breaks = seq(-2, plotLimit, histBreaksSize), plot = 0)), n = n()) %>% pull(a) %>%
                        #map(~Hist2tibble(.x)) %>% reduce(bind_rows)


#aggregating 
#SimAnnWetDays = SimAnnWetDays %>% group_by(breaks) %>% summarise(simulated = mean(density), min = quantile(density, 0.05), median = quantile(density, 0.5), max = quantile(density, 0.95)) %>%
                #left_join(IMSAnnWetDays, by = "breaks")

p3 = ggplot(IMSAnnWetDaysE %>% left_join(IMSAnnWetDaysS, by = "breaks"), aes(x = breaks)) +
                                ylab("pdf") +
                                geom_line(aes(y = elat, color = "Elat")) +
                                geom_line(aes(y = sdom, color = "Sedom")) +
                                xlab("# Annual wet days") +
                                     scale_y_continuous(expand = c(0, 0.0)) + scale_x_continuous(breaks = seq(0, plotLimit, 5), expand = c(0, 0));








#PET density----
plotLimit = ceiling(max(SynthRain$PET));
histBreaksSize = 0.1;

#histogram for observed
DailyPETE = elat[[2]] %>% filter(!is.na(pen)) %>% pull(pen) %>%
                                                        hist(breaks = seq(0, plotLimit, histBreaksSize), plot = 0) %>% Hist2tibble() %>% rename(elat = density)
DailyPETS = sedom[[2]] %>% filter(!is.na(pen)) %>% pull(pen) %>%
                                                        hist(breaks = seq(0, plotLimit, histBreaksSize), plot = 0) %>% Hist2tibble() %>% rename(sdom = density)

#SimPet = SynthRain %>% group_by(PETSeries) %>%
                       #summarise(a = list(hist(PET, breaks = seq(0, plotLimit, histBreaksSize), plot = 0)), n = n()) %>% pull(a) %>%
                             #map(~Hist2tibble(.x)) %>% reduce(bind_rows)
#SimPetDensity = SimPet %>% group_by(breaks) %>% summarise(simulated = mean(density), min = quantile(density, 0.05), median = quantile(density, 0.5), max = quantile(density, 0.95)) %>%
                #left_join(DailyPET, by = "breaks") %>% left_join(DailyPETWet, by = "breaks") %>% left_join(DailyPETDAD, by = "breaks") %>% left_join(DailyPETDAW, by = "breaks")

p5 = ggplot(DailyPETE %>% left_join(DailyPETS, by = "breaks"), aes(x = breaks)) +
                                ylab("pdf") +
                                 geom_line(aes(y = elat, color = "Elat")) +
                                 geom_line(aes(y = sdom, color = "Sedom")) +
                                 xlab("PET [mm / day]") +
                                 scale_y_continuous(expand = c(0, 0.0)) + scale_x_continuous(limits = c(0, plotLimit), expand = c(0, 0));


