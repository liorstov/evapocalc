#print gypsum path
draw = list(test1,test2,test3,test4) %>% map_dfc(~tibble(gyp = .x$YearGyp, gyp_depth = .x$gypDepth,module = factor(.x$module)))  %>% rowid_to_column(var = "year")#%>% mutate_all(~MovingAvarage(.x, 100, 100))
draw %>% gather("profile", "totalGypsum", - year,-module) %>% mutate(categ = ifelse(str_detect(profile, "depth"), "depth", "concentration")) %>% ggplot(aes(year, totalGypsum, color = (module))) + coord_cartesian(ylim = c()) + geom_point() + scale_y_continuous(name = "total gypsum [meq]") + xlab("duration [yr]") + guides(color = guide_legend(override.aes = list(size = 8))) + labs(color = "with FC module") + facet_wrap(strip.position = NULL, ncol = 2, categ ~ .,scales = "free_y") + scale_color_brewer(type = "qual", palette = 2)

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


draw = list(test1, test2, test3, test4) %>% map_dfr(~tibble(year = 1:length(.x$YearGyp), gyp = .x$YearGyp, gyp_depth = .x$gypDepth, module = (.x$module)))
a = draw %>% ggplot(aes(year, gyp, color = (module))) + coord_cartesian(ylim = c()) + geom_point() + scale_y_continuous(name = "total gypsum\n[meq/100 gr soil]") + xlab("duration [yr]") + guides(color = guide_legend(override.aes = list(size = 8))) + labs(color = "with FC module") + scale_color_brewer(type = "qual", palette = 2) + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
b = draw %>% ggplot(aes(year, gyp_depth, color = (module)),alpha = 0.8) + coord_cartesian(ylim = c()) + geom_point() + scale_y_reverse(name = "gypsum\ndepth [cm]") + xlab("duration [yr]") + guides(color = guide_legend(override.aes = list(size = 8))) + labs(color = "with FC module") + scale_color_brewer(type = "qual", palette = 2) + theme(legend.position =  "none")
grid.arrange(ggarrange(a, b)) 

draw = list(test1, test2, test3, test4) %>% map_dfr(~tibble(.x$WD[, 1:2], module = (.x$module))) %>% mutate(maxWD = MovingAvarage(maxWD, 500, 500))
draw  %>% ggplot(aes(year, maxWD, color = (module))) + geom_point() +coord_cartesian(xlim =  c())+ scale_y_reverse(name = "annual maximad WD [cm]") + xlab("duration [yr]") + guides(color = guide_legend(override.aes = list(size = 8))) + labs(color = "")  + scale_color_brewer(type = "qual", palette = 2)
