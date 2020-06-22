setwd('/home/username_here/PersonalizedSchedules')
rm(list = ls())

library(splines)
if(!require(JMbayes)){
  install.packages("JMbayes")
}

if(!require(MASS)){
  install.packages("MASS")
}

#Load a model that is fitted to the original PRIAS dataset.
#using this model we will generate data that tries to mimic PRIAS dataset
load("data/fitted_model_original_prias_dataset.Rdata")

source("src/utility_files/VariousSchedules.R")
source("src/utility_files/SyntheticDataSettingsAndFunctions.R")

args = commandArgs(trailingOnly=TRUE)
#number of patients for synthetic data
nSub = as.numeric(args[1])

#seed selected based on the current year 2020
set.seed(2020)

#This call generates PSA data, DRE data and event times
synthetic_data = generateSyntheticData(nSub = nSub)

#the synthetic_data object returns longitudinal data in simDs sub-object 
#and returns time to event data in simDs.id sub-object
longitudinal_data = synthetic_data$simDs
longitudinal_data$progression_time = rep(synthetic_data$simDs.id$progression_time, each=synthetic_data$timesPerSubject)
longitudinal_data$progressed = rep(synthetic_data$simDs.id$progressed, each=synthetic_data$timesPerSubject)

#Some patients will have loss to follow-up
#We assume that loss to follow-up happens only after year one of follow-up, and up to a maximum time of 15 years
#As long as loss to follow-up is non-informative censoring, changing censoring distribution does not affect model parameter estimation
loss_to_follow_up_max_time = 15

#Loss to folllow-up distribution. Even if you change, please note loss to follow-up before year one is very rare.
loss_to_follow_up_times = runif(n = nSub, min = 1, max = loss_to_follow_up_max_time)
longitudinal_data$loss_to_follow_up_time = rep(loss_to_follow_up_times, each=synthetic_data$timesPerSubject)
longitudinal_data = longitudinal_data[longitudinal_data$year_visit <= longitudinal_data$loss_to_follow_up_time,]

#Now we generate interval censored event times, using PRIAS based biopsy schedule
longitudinal_data = do.call('rbind', lapply(split(longitudinal_data, f=longitudinal_data$P_ID), FUN = function(pat_df){
  latest_survival_time = 0
  earliest_failure_time = Inf
  true_progression_time = pat_df$progression_time[1]
  progressed = pat_df$progressed[1]
  
  for(i in 1:nrow(pat_df)){
    current_visit_time = pat_df$year_visit[i]
    if(ifPRIASBiopsy(pat_df[1:i,], current_visit_time, latest_survival_time)){
      if(current_visit_time >= true_progression_time & progressed==T){
        earliest_failure_time = current_visit_time
        break
      }else{
        latest_survival_time = current_visit_time
      }
    }
  }
  
  pat_df$latest_survival_time = latest_survival_time
  pat_df$earliest_failure_time = earliest_failure_time
  
  return(pat_df)
}))

longitudinal_data = longitudinal_data[,c("P_ID", "age", "year_visit", 
                                         "log2psaplus1", "palpable_dre",
                                         "latest_survival_time", "earliest_failure_time")]

rm(list=setdiff(ls(), "longitudinal_data"))

save(longitudinal_data, file='data/synthetic_data.Rdata')