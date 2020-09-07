#print gypsum path
test %>% map_dfc(~tibble(gyp = .x$YearGyp, maxGyp = .x$YearMaxGyp * 8, gyp_depth = (.x$gypDepth + 0.5) * 40)) %>% mutate_all(~MovingAvarage(.x, 100, 100)) %>% rowid_to_column(var = "year") %>% gather("profile", "totalGypsum", - year) %>% ggplot(aes(year, totalGypsum, color = profile)) + coord_cartesian(ylim = c()) + geom_point() + scale_y_continuous(name = "total gypsum [meq]", sec.axis = sec_axis(~. / 8, name = "Gypsic Horizon [meq] / horizon depth[cm]")) + xlab("duration [yr]") + guides(color = guide_legend(override.aes = list(size = 8)))
test1 %>% list %>% map_dfc(~tibble(runoff = 100 * (1 - .x$rainStat$annual / .x$rainStat$origAnnual))) %>% mutate_all(~MovingAvarage(.x, 100, 100)) %>% rowid_to_column(var = "year") %>% ggplot(aes(year, runoff)) + coord_cartesian(ylim = c()) + geom_point() + scale_y_continuous(name = "runoff [%]") + xlab("duration [yr]")

bla = results$T1.10[c(1)] %>% map_dfc(~tibble(.x$rainStat, gyp = .x$YearGyp, maxgyp = .x$YearMaxGyp, so4 = .x$YearSulfate, ca = .x$YearCa)) %>% mutate(gypagg = gyp - lag(gyp), caagg = ca - lag(ca), so4agg = so4 - lag(so4))
bla %>% dplyr::select(gyp, gyp1, year, PET, PET1) %>% mutate_all(~MovingAvarage(.x, 100)) %>% gather("profile", "value", - year) %>%
    ggplot(aes(year, value, color = profile)) + geom_point() + labs(y = "[mm] / [meq]")

x = bla %>% filter(gypagg > 0) %>% pull(gypagg)
qnt <- quantile(x, probs = c(.25, .75), na.rm = T)
caps <- quantile(x, probs = c(.05, .95), na.rm = T)
H <- 1.5 * IQR(x, na.rm = T)
bla %>% filter(gypagg > (qnt[1] - H), gypagg < (qnt[2] + H)) %>% ggplot(aes(annual, gypagg)) + geom_point() + coord_cartesian(ylim = c()) + scale_color_gradientn(colours = rainbow(5, rev = T)) + labs(y = "annual gypsum percipitation [meq/100 gr soil]", x = "annual rain [mm]")
bla %>% filter(gypagg < 0) %>% ggplot(aes(annual, gypagg)) + geom_point() + coord_cartesian(ylim = c()) + scale_color_gradientn(colours = rainbow(5, rev = T)) + labs(y = "annual gypsum percipitation [meq/100 gr soil]", x = "annual rain [mm]")
bla %>% ggplot(aes(annual, gypagg, color = n)) + geom_point() + coord_cartesian(ylim = c()) + scale_color_gradientn(colours = rainbow(5, rev = T))

 colnames(bla) = c("gyp", "maxGyp", "gyp_depth", "gyp_NoRunoff", "maxGyp_NoRunoff", "gyp_depth_NoRunoff")
