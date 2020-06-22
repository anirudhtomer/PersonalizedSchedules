#Only works for one patient at a time
#object is the joint model object 
#new data is patient data
#last_test_time is time of last negative test on which progression was not detected
#planned_test_schedule is the schedule that has to be evaluated
#seed: random seed 
#M number of MCMC samples
#cache_size controls the accuracy of the results. higher cache size may give more accurate results, but will also take more time
getBurdenAndBenefit <- function (object, newdata, idVar = "id", last_test_time = NULL,
                                 planned_test_schedule, seed = 1L, M = 400L, cache_size=1000,...) {
  
  if(length(unique(newdata[[idVar]])) > 1){
    stop("Expected time delay and number of tests can be calculated for only one patient at a time.\n")
  }
  timeVar <- object$model_info$timeVar
  
  if(is.null(last_test_time)){
    last_test_time <- max(newdata[[timeVar]], na.rm = T)
  }
  
  cur_visit_time = max(newdata[[timeVar]], na.rm = T)
  
  planned_test_schedule = planned_test_schedule[planned_test_schedule >= cur_visit_time]
  if(length(planned_test_schedule)==0){
    stop("Expected time delay and number of tests cannot be calculated for an empty schedule or a schedule with maximum visit time less than current visit time.\n")
  }
  
  #Step 1: Create a CACHE of predicted survival probabilities to speed up operation
  SURV_CACHE_TIMES <- seq(from = last_test_time, to = max(planned_test_schedule), length.out = cache_size)
  SURV_CACHE_TIMES <- unique(sort(c(planned_test_schedule, SURV_CACHE_TIMES), decreasing = F))
  
  surv_pred <- survfitJM(object=object, newdata=newdata, idVar=idVar, survTimes = SURV_CACHE_TIMES[-1],
                         last.time = last_test_time, seed=seed, M=M, ...)
  
  SURV_CACHE_FULL <- rbind(rep(1, M), 
                           do.call('cbind', lapply(surv_pred$full.results, "[[", 1)))
  
  SURV_CACHE_FULL_rescaled = 1 - t(apply(1-SURV_CACHE_FULL, 1, FUN = function(x){
    x / (1-SURV_CACHE_FULL[nrow(SURV_CACHE_FULL),])
  })) 
  
  wk = JMbayes:::gaussKronrod()$wk
  sk = JMbayes:::gaussKronrod()$sk
  getGaussianQuadWeightsPoints = function(lower_upper_limit){
    p1 = 0.5 * sum(lower_upper_limit)
    p2 = 0.5 * diff(lower_upper_limit)
    
    return(list(weights = p2*wk, points = p2 * sk + p1))
  }
  
  ###############################################
  # Consequences
  ###############################################
  consequences = function(test_times){
    
    #Now we create intervals in which to integrate the conditional surv prob
    if(length(test_times)==1){
      test_intervals = list(c(last_test_time, test_times))
    }else{
      test_intervals = c(last_test_time,
                         rep(test_times[-length(test_times)],each=2),
                         test_times[length(test_times)])
      test_intervals = split(test_intervals, rep(1:(length(test_intervals)/2), each=2))
    }
    
    interval_res = vector("list", length(test_intervals))
    exp_num_tests = rep(0, M)
    for(j in 1:length(test_intervals)){
      wt_points=getGaussianQuadWeightsPoints(test_intervals[[j]])
      lower_limit = test_intervals[[j]][1]
      upper_limit = test_intervals[[j]][2]
      
      lower_limit_nearest_index = which(lower_limit==SURV_CACHE_TIMES)
      upper_limit_nearest_index = which(upper_limit==SURV_CACHE_TIMES)
      interval_res[[j]]$cum_risk_interval = SURV_CACHE_FULL_rescaled[lower_limit_nearest_index,] - 
        SURV_CACHE_FULL_rescaled[upper_limit_nearest_index,]
      
      cond_expected_fail_time = sapply(1:length(wt_points$points), function(i){
        cum_surv_at_points = SURV_CACHE_FULL_rescaled[which.min(abs(wt_points$points[i]-SURV_CACHE_TIMES)),] - SURV_CACHE_FULL_rescaled[upper_limit_nearest_index,]
        scaled_cum_surv_at_points = cum_surv_at_points/interval_res[[j]]$cum_risk_interval
        
        return(wt_points$weights[i] * scaled_cum_surv_at_points)
      })
      
      interval_res[[j]]$delay = upper_limit - (lower_limit + apply(cond_expected_fail_time, 1, sum))
      
      exp_num_tests = exp_num_tests + j * interval_res[[j]]$cum_risk_interval
    }
    exp_num_tests = exp_num_tests + j * SURV_CACHE_FULL_rescaled[upper_limit_nearest_index,]
    exp_num_tests = mean(exp_num_tests, na.rm = T)
    
    expected_detection_delay = mean(apply(sapply(interval_res, FUN = function(x){
      x$delay * x$cum_risk_interval
    }),1, sum, na.rm=T),na.rm=T)
    
    #proposed test times always has a final test at year 10
    return(list(expected_detection_delay = expected_detection_delay,
                expected_num_tests = exp_num_tests,
                planned_test_schedule = test_times))
  }
  
  test_schedule_conseq = consequences(planned_test_schedule)
  test_schedule_conseq$euclidean_distance = sqrt(test_schedule_conseq$expected_detection_delay^2 + (test_schedule_conseq$expected_num_tests-1)^2)
  
  return(test_schedule_conseq)
}

