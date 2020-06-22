fixed_psaFormula = ~ 1 + I(age - 65) +  I((age - 65)^2) + ns(year_visit, knots=c(0.75, 2.12), Boundary.knots=c(0, 6.4))
random_psaFormula = ~ 1 + ns(year_visit, knots=c(0.75, 2.12), Boundary.knots=c(0, 6.4))
fixed_random_psaSlopeFormula = ~ 0 + dns(year_visit, knots=c(0.75, 2.12), Boundary.knots=c(0, 6.4))

fixed_dreFormula = ~ 1 + I(age - 65) +  I((age - 65)^2) + year_visit
random_dreFormula = ~ 1 + year_visit

survivalFormula = ~ 0 + I(age - 65) + I((age - 65)^2)

bNames = c("b_Int_DRE", "b_Slope_DRE", "b_Int_PSA", "b_Slope1_PSA", "b_Slope2_PSA", "b_Slope3_PSA")

#t-distibution (3 degress of freedom) error terms
psaErrorDist = "t3"
#data generated only until year ten
MAX_FAIL_TIME = 10

IND_RANDOM_DRE = 1:2
IND_RANDOM_PSA = 3:6
IND_RANDOM_PSA_SLOPE = IND_RANDOM_PSA[-1]
IND_FIXED_PSA_SLOPE = 4:6

forms_simDs = list("palpable_dre"= "value",
                   "log2psaplus1" = "value",
                   "log2psaplus1" = list(fixed =fixed_random_psaSlopeFormula,
                                         random=fixed_random_psaSlopeFormula,
                                         indFixed = 4:6, indRandom=2:4, name = "slope"))

generatePSATimeBySchedule = function(){
  months = c(seq(0, 24, 3), seq(30, MAX_FAIL_TIME*12, 6))
  
  return(months/12)
}

generateDRETimeBySchedule = function(){
  months = c(seq(0, MAX_FAIL_TIME*12,6))
  return(months/12)
}

generateTDistError = function(n){
  sigma_psa = fitted_model_original_prias_dataset$statistics$postMeans$sigma2
  JMbayes:::rgt(n = n, mu = 0, sigma = sigma_psa, df = 3)
}

generateNormDistError = function(n){
  sigma_psa = fitted_model_original_prias_dataset$statistics$postMeans$sigma2
  rnorm(n=n, mean = 0, sd = sigma_psa)
}

generateTruePSAProfile = function(age, year_visit, randomEff_psa){
  betas_psa = fitted_model_original_prias_dataset$statistics$postMeans$betas2
  df = data.frame(age, year_visit)
  model.matrix(fixed_psaFormula, df) %*% betas_psa + model.matrix(random_psaFormula, df) %*% as.numeric(randomEff_psa)
}

generateTruePSASlope = function(year_visit, randomEff_psa_slope){
  betas_psa_time = fitted_model_original_prias_dataset$statistics$postMeans$betas2[IND_FIXED_PSA_SLOPE]
  df = data.frame(year_visit)
  model.matrix(fixed_random_psaSlopeFormula, df) %*% betas_psa_time + model.matrix(fixed_random_psaSlopeFormula, df) %*% as.numeric(randomEff_psa_slope)
}

generateTrueDRELogOdds = function(age, year_visit, randomEff_dre){
  betas_dre = fitted_model_original_prias_dataset$statistics$postMeans$betas1
  df = data.frame(age, year_visit)
  model.matrix(fixed_dreFormula, df) %*% betas_dre + model.matrix(random_dreFormula, df) %*% as.numeric(randomEff_dre)
}

generatePalpableDRE = function(logOdds){
  rbinom(length(logOdds), 1, plogis(logOdds))
}

generateSyntheticData = function(nSub, psaErrorDist = "t3"){
  
  print("Please be patient: it takes 1 second for 5 patients on a 4th generation Intel core-i7 processor")
  
  survivalFunc <- function (t, patientId) {
    exp(-integrate(hazardFunc, lower = 0, upper = t, patientId)$value)
  }
  
  invSurvival <- function (t, u, patientId) {
    log(u) + integrate(hazardFunc, lower = 0, upper = t, patientId)$value
  }
  
  hazardFunc = function (year_visit, patientId) {
    alphas = fitted_model_original_prias_dataset$statistics$postMeans$alphas
    
    b_subject = simDs.id[patientId, bNames]
    age = simDs.id$age[patientId]  
    wGamma = simDs.id$wGamma[patientId]
    
    baselinehazard = exp(splineDesign(fitted_model_original_prias_dataset$control$knots, year_visit, 
                                      ord = fitted_model_original_prias_dataset$control$ordSpline, outer.ok = T) %*% fitted_model_original_prias_dataset$statistics$postMeans$Bs_gammas)
    
    truePSA = generateTruePSAProfile(age, year_visit, b_subject[IND_RANDOM_PSA])
    truePSASlope = generateTruePSASlope(year_visit, b_subject[IND_FIXED_PSA_SLOPE])
    trueDRELogOdds = generateTrueDRELogOdds(age, year_visit, b_subject[IND_RANDOM_DRE])
    
    y_Alpha = cbind(trueDRELogOdds, truePSA, truePSASlope) %*% alphas
    
    baselinehazard * exp(wGamma + y_Alpha)
  }
  
  pSurvTime = function(survProb, patientId){
    Low = 1e-05
    Up <- MAX_FAIL_TIME
    Root <- try(uniroot(invSurvival, interval = c(Low, Up), 
                        u = survProb, patientId)$root, silent = TRUE)
    if(inherits(Root, "try-error")){
      #print(Root)
      return(NA)
    }else{
      return(Root)
    }
  }
  
  longTimes <- c(replicate(nSub, generatePSATimeBySchedule(), simplify = T)) 
  timesPerSubject = length(longTimes) / nSub
  
  ###########################
  # Generate the longitudinal dataset
  ###########################
  subId <- 1:nSub
  
  D = fitted_model_original_prias_dataset$statistics$postMeans$D
  b <- mvrnorm(nSub, mu = rep(0, nrow(D)), D)
  colnames(b) = bNames
  
  age <- rnorm(n = nSub, mean=65.67, sd=6.821431)
  simDs.id = data.frame(P_ID = subId, age = age, progression_time = NA)
  
  gammas = fitted_model_original_prias_dataset$statistics$postMeans$gammas
  simDs.id$wGamma <- as.numeric(model.matrix(survivalFormula, data = simDs.id) %*% gammas)
  
  simDs.id = cbind(simDs.id, b)
  rm(b)
  
  #Longitudinal dataset
  simDs = data.frame(P_ID = rep(subId, each=timesPerSubject), 
                     age = rep(age, each=timesPerSubject),
                     year_visit = longTimes,
                     visitNumber = rep(1:timesPerSubject, nSub),
                     log2psaplus1 = NA, palpable_dre=NA)
  
  simDs$trueLogOddsPalpableDRE = unlist(by(data=simDs, simDs$P_ID, function(data){
    pid = data$P_ID[1]
    randomEff = simDs.id[pid, bNames[IND_RANDOM_DRE]]
    generateTrueDRELogOdds(age = data$age, year_visit = data$year_visit, 
                           randomEff_dre = randomEff)
  }))
  
  simDs$palpable_dre = generatePalpableDRE(logOdds = simDs$trueLogOddsPalpableDRE)
  simDs$palpable_dre[!(simDs$year_visit %in% generateDRETimeBySchedule())] = NA
  
  simDs$trueLog2psaplus1Profile = unlist(by(data=simDs, simDs$P_ID, function(data){
    pid = data$P_ID[1]
    randomEff = simDs.id[pid,bNames[IND_RANDOM_PSA]]
    generateTruePSAProfile(age = data$age, year_visit = data$year_visit, randomEff_psa = randomEff)
  }))
  
  if(psaErrorDist=="normal"){
    simDs$log2psaplus1 = simDs$trueLog2psaplus1Profile + generateNormDistError(nrow(simDs))
  }else{
    simDs$log2psaplus1 = simDs$trueLog2psaplus1Profile + generateTDistError(nrow(simDs))
  }
  
  ###################
  #Now generate event times
  ###################
  print("Done generating longitudinal data. Now trying to generate event times.")
  
  simDs.id$survProbs <- runif(nSub)
  
  simDs.id$progression_time = sapply(1:nSub, FUN = function(i){
    pSurvTime(simDs.id$survProbs[i], i)
  })
  
  percentageCensored = sum(is.na(simDs.id$progression_time))/nSub
  print(paste0("Event rate= ", 100-percentageCensored*100, "%"))
  
  simDs.id$progressed = ifelse(is.na(simDs.id$progression_time), 0, 1)
  simDs.id$studyend_censored = simDs.id$progressed
  simDs.id$progression_time[simDs.id$progressed==0] = MAX_FAIL_TIME
  
  return(list(simDs=simDs, simDs.id=simDs.id, timesPerSubject=timesPerSubject))
}

fitJointModelOnNewData = function(simDs, simDs.id, timesPerSubject, nSubTraining, nSubTest=0, 
                                  censStartTime=MAX_FAIL_TIME, censEndTime=MAX_FAIL_TIME){
  
  #Divide into training and test
  trainingDs.id = simDs.id[1:nSubTraining, ]
  testDs.id = simDs.id[(nSubTraining+1):(nSubTraining + nSubTest),]
  trainingDs = simDs[simDs$P_ID %in% trainingDs.id$P_ID, ]
  testDs = simDs[simDs$P_ID %in% testDs.id$P_ID, ]
  
  #Dropout censoring for training patients
  Ctimes<-runif(nSubTraining, censStartTime, censEndTime)
  
  trainingDs.id$progressed = trainingDs.id$progressed==1 & trainingDs.id$progression_time <= Ctimes
  trainingDs.id$progression_time = pmin(trainingDs.id$progression_time, Ctimes)
  trainingDs.id$Ctime = Ctimes
  
  trainingDs$Ctime = rep(trainingDs.id$Ctime, each=timesPerSubject)
  trainingDs$progression_time = rep(trainingDs.id$progression_time, each=timesPerSubject)
  trainingDs$progressed = rep(trainingDs.id$progressed, each=timesPerSubject)
  testDs$progression_time = rep(testDs.id$progression_time, each=timesPerSubject)
  
  # drop the longitudinal measurementsthat were taken after the observed event time for each subject.
  trainingDs = trainingDs[trainingDs$year_visit <= trainingDs$progression_time, ]
  
  print("Starting to fit long model with mvglmer")
  
  startTime_mvglmer_simDs = Sys.time()
  mvglmer_dre_psa_simDs = mvglmer(list(palpable_dre~I(age - 65) + I((age - 65)^2) + year_visit  +
                                         (year_visit|P_ID),
                                       
                                       log2psaplus1 ~ I(age - 65) + I((age - 65)^2) +
                                         ns(year_visit, knots=c(0.75, 2.12), Boundary.knots=c(0, 6.4)) +
                                         (ns(year_visit, knots=c(0.75, 2.12), Boundary.knots=c(0, 6.4))|P_ID)),
                                  
                                  data=trainingDs, families = list(binomial, gaussian), engine = "STAN",
                                  control = list(n.iter=1000, n.processors=max_cores))
  endTime_mvglmer_simDs = Sys.time()
  mvglmer_fitting_time = endTime_mvglmer_simDs - startTime_mvglmer_simDs
  
  print("Done fitting long model with mvglmer")
  print(mvglmer_fitting_time)
  
  survModel_simDs = coxph(Surv(progression_time, progressed) ~ I(age - 65) + I((age - 65)^2),
                          data = trainingDs.id, x=T, model=T)
  
  print("Starting to fit joint model with mvJointModelBayes")
  startTime_mvJoint_simDs = Sys.time()
  mvJoint_dre_psa_simDs = mvJointModelBayes(mvglmer_dre_psa_simDs, survModel_simDs,
                                            timeVar = "year_visit", Formulas = forms_simDs,
                                            control=list(n_cores=max_cores))
  endTime_mvJoint_simDs = Sys.time()
  mvjoint_fitting_time = endTime_mvJoint_simDs - startTime_mvJoint_simDs
  
  print("Done fitting joint model with mvJointModelBayes")
  print(mvjoint_fitting_time)
  
  out = list("trainingData"=list(trainingDs=trainingDs, trainingDs.id=trainingDs.id),
             "testData"=list(testDs=testDs, testDs.id=testDs.id),
             "censStartTime"=censStartTime, "censEndTime"=censEndTime, 
             "survModel_simDs"=survModel_simDs, "mvglmer_dre_psa_simDs"=mvglmer_dre_psa_simDs, 
             "mvglmer_fitting_time"=mvglmer_fitting_time,
             "mvJoint_dre_psa_simDs"=mvJoint_dre_psa_simDs, 
             "mvjoint_fitting_time"=mvjoint_fitting_time)
  
  return(out)
}
