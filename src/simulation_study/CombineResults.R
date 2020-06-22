setwd('/home/username_here/PersonalizedSchedules')

args = commandArgs(trailingOnly=TRUE)
seed_range = seq(as.numeric(args[1]), as.numeric(args[2]))

simstudy_results = do.call('rbind', lapply(seed_range, function(seed){
  files = list.files("data/simulation_study/schedule_results/",
                     pattern = as.character(seed), full.names = T)
  
  biopsyDfs = lapply(files, FUN = function(file){
    load(file)
    return(biopsyDf)
  })
  
  biopsyDf_summary_list = lapply(biopsyDfs, FUN = function(patient_df){
    schedules = unique(patient_df$schedule)
    
    nb = sapply(schedules, FUN = function(schedule){
      nrow(patient_df[patient_df$schedule==schedule,])
    })
    
    delay = sapply(schedules, FUN = function(schedule){
      tail(patient_df$biopsy_time[patient_df$schedule==schedule],1) - patient_df$progression_time[1]
    })
    
    ret_df = patient_df[rep(1, length(schedules)), c("seed", "P_ID", "progression_time", "progressed")]
    ret_df$schedule = schedules
    ret_df$nb = nb
    ret_df$delay = delay
    return(ret_df)
  })
  
  biopsyDf_summary = do.call('rbind', biopsyDf_summary_list)
  return(biopsyDf_summary)
}))

save(simstudy_results, file="data/simulation_study/simstudy_results.Rdata")