setwd('/home/username_here/PersonalizedSchedules')

args = commandArgs(trailingOnly=TRUE)

seed = as.numeric(args[1])
testId = as.numeric(args[2])

library(splines)
if(!require(JMbayes)){
  install.packages("JMbayes")
}

if(!require(MASS)){
  install.packages("MASS")
}

if(!require(survival)){
  install.packages("survival")
}

load(paste0("data/simulation_study/fitted_models/jointModelData_seed_", seed, ".Rdata"))
source("src/utility_files/VariousSchedules.R")

ANNUAL = "Annual"
PRIAS = "PRIAS"
FixedRisk_05 = "FixedRisk_5%"
FixedRisk_10 = "FixedRisk_10%"
AutoRisk_Delay_75 = "AutoRisk_Delay_0.75"
AutoRisk_NoDelay = "AutoRisk_NoDelayLimit"

#M is the number of MCMC samples
M=500
#We conduct biopsies until year ten of follow-up
MAX_FAIL_TIME = 10

#First biopsy is decided at visit number 5, that is, year one of follow-up
FIRST_DECISION_VISIT_NR = 5

#minimum gap between biopsies is one year (recommended by almost all prostate cancer surveillance programs)
MIN_BIOPSY_GAP = 1

testDs = jointModelData$testData$testDs[jointModelData$testData$testDs$P_ID == testId,]
testDs.id = jointModelData$testData$testDs.id[jointModelData$testData$testDs.id$P_ID == testId,]
progression_time = testDs.id$progression_time

set.seed(seed)

print(paste("**** Patient ID: ", testId, 
            " with progression time: ", progression_time,
            " from simulation with seed: ", seed))

#First we run annual and prias biopsy schedules
if(progression_time<=1){
  annual_biopsies = prias_biopsies = 1
}else{
  annual_biopsies = seq(1, ceiling(progression_time), by = 1)
  #PRIAS schedule
  prias_biopsies = c(1)
  for(row_num in FIRST_DECISION_VISIT_NR:nrow(testDs)){
    patient_df = testDs[1:row_num,]
    cur_visit_time = testDs$year_visit[row_num]
    
    if(ifPRIASBiopsy(patient_df, cur_visit_time, max(prias_biopsies,0))){
      prias_biopsies = c(prias_biopsies, cur_visit_time)
    }
    
    if(max(prias_biopsies, 0) >= progression_time){
      break
    }
  }
  #check if progression could not be detected within 10 years despite being progressed
  if(max(prias_biopsies, 0) < progression_time){
    prias_biopsies = c(prias_biopsies, MAX_FAIL_TIME)
  }
}

print('done annual and prias schedules')

#Next we do risk based biopsies based on fixed thresholds
fixed_risk_biopsies = lapply(c(0.05, 0.1), function(threshold){
  source("src/utility_files/prediction_psa_dre.R")
  risk_biopsies = c()
  for(row_num in FIRST_DECISION_VISIT_NR:nrow(testDs)){
    patient_df = testDs[1:row_num,]
    cur_visit_time = testDs$year_visit[row_num]
    
    if(ifFixedRiskBasedBiopsy(object = jointModelData$mvJoint_dre_psa_simDs,
                              patient_data = patient_df, 
                              cur_visit_time = cur_visit_time, 
                              last_biopsy_time = max(risk_biopsies, 0),
                              min_biopsy_gap = MIN_BIOPSY_GAP,
                              threshold = threshold,
                              M = M)){
      risk_biopsies = c(risk_biopsies, cur_visit_time)
    }
    
    if(max(risk_biopsies, 0) >= progression_time){
      break
    }
  }
  if(max(risk_biopsies, 0) < progression_time){
    risk_biopsies = c(risk_biopsies, MAX_FAIL_TIME)
  }
  
  print(paste0('done fixed risk biopsy with threshold: ', threshold*100, "%"))
  return(risk_biopsies)
})

#Finally we conduct biopsies using a dynamically chosen risk threshold 
#based on the utility function proposed in our paper
auto_risk_biopsies = lapply(c(0.75, Inf), function(delay_limit){
  source("src/utility_files/prediction_psa_dre_randEff_reuse.R")
  automatic_kappa_biopsies = c()
  for(row_num in FIRST_DECISION_VISIT_NR:nrow(testDs)){
    patient_df = testDs[1:row_num,]
    cur_visit_time = testDs$year_visit[row_num]
    
    if(ifAutomaticRiskBasedBiopsy(object = jointModelData$mvJoint_dre_psa_simDs,
                                  patient_data = patient_df, 
                                  cur_visit_time = cur_visit_time, 
                                  last_biopsy_time = max(automatic_kappa_biopsies, 0),
                                  min_biopsy_gap = MIN_BIOPSY_GAP, delay_limit = delay_limit,
                                  M = M, horizon = MAX_FAIL_TIME, use_exact_delay = T)){
      automatic_kappa_biopsies = c(automatic_kappa_biopsies, cur_visit_time)
    }
    
    if(max(automatic_kappa_biopsies, 0) >= progression_time){
      break
    }
  }
  if(max(automatic_kappa_biopsies, 0) < progression_time){
    automatic_kappa_biopsies = c(automatic_kappa_biopsies, MAX_FAIL_TIME)
  }
  print(paste0('done automatic risk biopsy with delay: ', delay_limit, " years"))
  return(automatic_kappa_biopsies)
})

scheduleNames = c(ANNUAL, PRIAS, 
                  FixedRisk_05, FixedRisk_10, 
                  AutoRisk_Delay_75, AutoRisk_NoDelay)

scheduleLengths = c(length(annual_biopsies), length(prias_biopsies),
                    sapply(fixed_risk_biopsies, length), sapply(auto_risk_biopsies, length))

biopsy_time = c(annual_biopsies, prias_biopsies,
                do.call('c', fixed_risk_biopsies), do.call('c', auto_risk_biopsies))

#Now we combine all of these into a data frame for this patient
biopsyDf = data.frame(seed = seed, P_ID = testId, 
                      progression_time = progression_time,
                      progressed = testDs.id$progressed,
                      schedule = rep(scheduleNames,scheduleLengths),
                      biopsy_time = biopsy_time)

print("Saving now")
patient_file_name = paste0("patient_seed_", seed, "_P_ID_", testId)
save(biopsyDf, file=paste0("data/simulation_study/schedule_results/", patient_file_name, ".Rdata"))