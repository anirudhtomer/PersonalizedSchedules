if(!require(ggplot2)){
  install.packages("ggplot2")
}

if(!require(ggpubr)){
  install.packages("ggpubr")
}

load(file = "data/simulation_study/simstudy_results.Rdata")

simstudy_results = droplevels(simstudy_results[!simstudy_results$schedule %in% c("FixedRisk_5%"),])

levels(simstudy_results$schedule)[3] = "\u03BA=10%"
levels(simstudy_results$schedule)[4] = "\u03BA*{v | E(D) \u2264 0.75}"
levels(simstudy_results$schedule)[5] = "\u03BA*(v)"

FONT_SIZE = 12
POINT_SIZE = 1.5
THEME_COLOR = 'dodgerblue4'
WARNING_COLOR = 'darkorange'

a = ggplot(simstudy_results[simstudy_results$progressed==1,], aes(x=schedule, y=nb)) +
  geom_boxplot(outlier.shape = NA) +
  stat_summary(fun.y=mean, geom="point", size=POINT_SIZE, color=WARNING_COLOR) +
  theme_bw() + 
  scale_y_continuous(breaks = c(1,4,7,10), 
                     limits = c(1,10), minor_breaks = 1:10) +
  theme(text=element_text(size=FONT_SIZE),
        title = element_text(size=FONT_SIZE-2),
        axis.title = element_text(size=FONT_SIZE)) +
  xlab("Schedule") + ylab("Number of biopsies") +
  coord_flip() +
  ggtitle("Progressing patients (50%)")

b = ggplot(simstudy_results[simstudy_results$progressed==1,], aes(x=schedule, y=delay)) +
  geom_boxplot(outlier.shape = NA)+
  stat_summary(fun.y=mean, geom="point", size=POINT_SIZE, color=WARNING_COLOR) +
  theme_bw() + 
  scale_y_continuous(breaks = 0:3, limits = c(0,3)) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        text=element_text(size=FONT_SIZE)) +
  xlab("Schedule") + ylab("Time delay in detecting\nprogression (years)")+
  coord_flip()

c = ggplot(simstudy_results[simstudy_results$progressed==0,], aes(x=schedule, y=nb)) +
  geom_boxplot(outlier.shape = NA) +
  stat_summary(fun.y=mean, geom="point", size=POINT_SIZE, color=WARNING_COLOR) +
  theme_bw() +
  scale_y_continuous(breaks = c(1,4,7,10), limits = c(1,10),
                     minor_breaks = 1:10) +
  theme(text=element_text(size=FONT_SIZE),
        title = element_text(size=FONT_SIZE-2),
        axis.title = element_text(size=FONT_SIZE)) +
  xlab("Schedule") + ylab("Number of biopsies") +
  coord_flip() +
  ggtitle("Non-progressing patients (50%)")

d = ggplot() + geom_text(aes(x=1, y=1.5), label="Time delay not available for\nnon-progressing patients")+
  theme_bw() + 
  scale_y_continuous(breaks = 0:3, limits = c(0,3)) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        text=element_text(size=FONT_SIZE), axis.ticks.y=element_blank()) +
  xlab("Schedule") + ylab("Time delay in detecting\nprogression (years)")+ 
  coord_flip()

upper_plot = ggarrange(a,b, nrow=1, ncol=2, align = "h", widths = c(1.5,1))
lower_plot = ggarrange(c,d, nrow=1, ncol=2, align = "h", widths = c(1.5,1))

final_plot = ggpubr::ggarrange(upper_plot, lower_plot, 
                  nrow=2, ncol=1, align = "v", labels = "AUTO")

print(final_plot)