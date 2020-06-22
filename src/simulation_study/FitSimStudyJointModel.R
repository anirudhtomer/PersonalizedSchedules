setwd('/home/username_here/PersonalizedSchedules')

args = commandArgs(trailingOnly=TRUE)

seed = as.numeric(args[1])
max_cores = as.numeric(args[2])

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

#Load a model that is fitted to the original PRIAS dataset.
#using this model we will generate data that tries to mimic PRIAS dataset
load("data/fitted_model_original_prias_dataset.Rdata")

source("src/utility_files/VariousSchedules.R")
source("src/utility_files/SyntheticDataSettingsAndFunctions.R")

print(paste("Using seed:", seed))
set.seed(seed)
newData = generateSyntheticData(nSub = 1000)
print('Generated data, now trying to fit the model')
jointModelData = fitJointModelOnNewData(simDs = newData$simDs, simDs.id=newData$simDs.id,
                                        timesPerSubject = newData$timesPerSubject,
                                        nSubTraining = 750, nSubTest = 250, 
                                        censStartTime = 1, censEndTime = 15)
jointModelData$seed = seed

print("********* Reducing model size before saving it ******")
jointModelData$mvglmer_dre_psa_simDs = NULL
jointModelData$mvglmer_fitting_time = NULL
jointModelData$mvjoint_fitting_time = NULL
jointModelData$survModel_simDs = NULL
jointModelData$mvJoint_dre_psa_simDs$mcmc$b = NULL
jointModelData$mvJoint_dre_psa_simDs$call = NULL
jointModelData$mvJoint_dre_psa_simDs$weights = NULL
jointModelData$mvJoint_dre_psa_simDs$mcmc_info = NULL
jointModelData$mvJoint_dre_psa_simDs$model_info = NULL
jointModelData$mvJoint_dre_psa_simDs$control = list(knots = jointModelData$mvJoint_dre_psa_simDs$control$knots,
                                                    ordSpline = jointModelData$mvJoint_dre_psa_simDs$control$ordSpline)

saveName = paste0("jointModelData_seed_",jointModelData$seed, ".Rdata")
save(jointModelData, file = paste0("data/simulation_study/fitted_models/", saveName))
