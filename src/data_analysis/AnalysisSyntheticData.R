setwd('/home/username_here/PersonalizedSchedules')
if(!require(survival)){
  install.packages("survival")
}

library(splines)
if(!require(JMbayes)){
  install.packages("JMbayes")
}

rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
load(args[1])

#number of parallel cores to use for fitting the joint model
max_cores = 2

#The joint model API requires the survival sub-model and longitudinal-sub model fitted separately

#First we fit a survival model for interval censored data. We use median centered age variable
survdata = longitudinal_data[!duplicated(longitudinal_data$P_ID),]
survdata$latest_survival_time[survdata$latest_survival_time==0] = 0.0001
survival_model = survreg(Surv(latest_survival_time, earliest_failure_time, type = "interval2") ~ 
                            I(age-65) +  I((age-65)^2), x = T, data = survdata, model = TRUE)


print("Started fitting long model with mvglmer")
#Then we fit a bivariate longitudinal data model. 
#palpable_dre is binary.
#log2psaplus1 is log_2(PSA + 1) transformed PSA values. we use splines to model evolution of PSA over time
startTime_mvglmer_simDs = Sys.time()
longitudinal_model = mvglmer(list(palpable_dre~I(age - 65) + I((age - 65)^2) + year_visit  +
                                       (year_visit|P_ID),
                                     
                                     log2psaplus1 ~ I(age - 65) + I((age - 65)^2) +
                                       ns(year_visit, knots=c(0.75, 2.12), Boundary.knots=c(0, 6.4)) +
                                       (ns(year_visit, knots=c(0.75, 2.12), Boundary.knots=c(0, 6.4))|P_ID)),
                                
                                data=longitudinal_data, families = list(binomial, gaussian), engine = "STAN",
                                control = list(n.iter=1000, n.processors=max_cores))
endTime_mvglmer_simDs = Sys.time()
mvglmer_fitting_time = endTime_mvglmer_simDs - startTime_mvglmer_simDs
print("Done fitting long model with mvglmer")
print(mvglmer_fitting_time)

#Finally the joint model is fitted 
print("Starting to fit joint model with mvJointModelBayes")
startTime_mvJoint_simDs = Sys.time()
fixed_random_psaSlopeFormula = ~ 0 + dns(year_visit, knots=c(0.75, 2.12), Boundary.knots=c(0, 6.4))

#Define the assocition between longitudinal and time-to-event (disease progression) outcomes.
forms_simDs = list("palpable_dre"= "value",
                   "log2psaplus1" = "value",
                   "log2psaplus1" = list(fixed =fixed_random_psaSlopeFormula,
                                         random=fixed_random_psaSlopeFormula,
                                         indFixed = 4:6, indRandom=2:4, name = "slope"))

fitted_jointmodel = mvJointModelBayes(longitudinal_model, survival_model,
                                          timeVar = "year_visit", Formulas = forms_simDs,
                                          control=list(n_cores=max_cores))
endTime_mvJoint_simDs = Sys.time()
mvjoint_fitting_time = endTime_mvJoint_simDs - startTime_mvJoint_simDs
print(mvjoint_fitting_time)

save(fitted_jointmodel, file = "data/fitted_jointmodel.Rdata")
