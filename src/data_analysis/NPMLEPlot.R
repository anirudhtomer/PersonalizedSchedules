if(!require(icenReg)){
  install.packages("icenReg")
}

load("data/synthetic_data.Rdata")

getSurvProbNPMLE = function(pred.time, baseline_est_int){
  index2 = tail(which(baseline_est_int[,2]<=pred.time),1)
  if(baseline_est_int[index2 + 1, 1]<=pred.time){
    pred.time.surv = baseline_est_int[index2 + 1,3]
  }else{
    pred.time.surv = baseline_est_int[index2, 3]
  }
  
  return(pred.time.surv)
}

npmle_prias = ic_sp(formula = Surv(time = latest_survival_time, time2 = earliest_failure_time,
                                   type = "interval2")~1, data = longitudinal_data[!duplicated(longitudinal_data$P_ID),])


baseline_est = getSCurves(npmle_prias)
baseline_est_int = baseline_est$Tbull_ints
baseline_est_int = cbind(baseline_est_int, baseline_est$S_curves$baseline)
baseline_est_int = rbind(c(0,0,1), baseline_est_int)

max_time = 10
timePoints = seq(0, max_time, length.out = 1000)
survProbs = sapply(timePoints, FUN = getSurvProbNPMLE, baseline_est_int=baseline_est_int)

plotdf = data.frame(timePoints=timePoints, riskProbs=1-survProbs)

FONT_SIZE=14
npmle_plot_all = ggplot() + 
  geom_line(aes(x=plotdf$timePoints, 
                y=plotdf$riskProbs))+
  scale_x_continuous(breaks=0:10, limits = c(0,10)) + 
  theme_bw() +
  theme(text = element_text(size=FONT_SIZE), 
        legend.position = "none",
        legend.text = element_text(size=FONT_SIZE-4),
        axis.line = element_line())+
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = paste0(seq(0, 1, 0.25)*100, "%"),
                     limits = c(0,1)) + 
  ylab("Cause-specific cumulative upgrading-risk (%)") +
  xlab("Follow-up time (years)")
print(npmle_plot_all)
