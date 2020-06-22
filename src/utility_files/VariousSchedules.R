#These schedules do not check if minimum gap is maintained between biopsies
ifFixedBiopsy = function(cur_visit_time, last_biopsy_time,
                         biopsy_frequency = 1){
  return(cur_visit_time - last_biopsy_time>=biopsy_frequency)
}

getFixedSchedule = function(cur_visit_time, biopsy_frequency=1, horizon=10){
  return(seq(cur_visit_time, horizon, by=biopsy_frequency))
}

#here we check if minimum gap is maintained
ifPRIASBiopsy = function(patient_data, cur_visit_time, last_biopsy_time){
  
  if(cur_visit_time - last_biopsy_time < 1){
    return(FALSE)
  }else{
    fixed_schedule = c(1, 4, 7, 10, 15)
    
    if(cur_visit_time %in% fixed_schedule){
      return(TRUE)
    }else{
      psa = pmax(0.1, 2^(patient_data$log2psaplus1)-1)
      log2psa = log(psa, base = 2)
      year_visit = patient_data$year_visit
      
      psa_dt = 1/(lm(log2psa~year_visit)$coefficients[2])
      
      return(psa_dt>=0 & psa_dt<=10)
    }
  }
}  

#assuming that observed longitudinal data is available up to current visit time
getPRIASSchedule = function(object, patient_data, cur_visit_time, 
                            last_biopsy_time,
                            M=500, horizon=10){
  
  PSA_CHECK_UP_TIME = c(seq(0, 2, 0.25), seq(2.5, 10, 0.5))
  visit_schedule = c(cur_visit_time, PSA_CHECK_UP_TIME[PSA_CHECK_UP_TIME > cur_visit_time & PSA_CHECK_UP_TIME <=horizon])
  future_log2psaplus1_matrix = getExpectedFutureOutcomes(object, patient_data, last_biopsy_time, 
                                                         long_predict_times = visit_schedule, 
                                                         psaDist = "Tdist", M = M, addRandomError = T)$predicted_psa
  future_log2psaplus1_matrix[1,] = rep(patient_data$log2psaplus1[nrow(patient_data)], M)
  patient_data = patient_data[-nrow(patient_data)]
  
  proposed_biopsy_times_list = vector("list", length=500)
  for(sim in 1:length(proposed_biopsy_times_list)){
    proposed_biopsy_times = c()
    temp_patient_data = patient_data
    previous_biopsy_time = last_biopsy_time
    
    for(i in 1:length(visit_schedule)){
      new_row = temp_patient_data[1,]
      new_row$year_visit = visit_schedule[i]
      new_row$log2psaplus1 = sample(x=future_log2psaplus1_matrix[i,], size=1)
      
      temp_patient_data = rbind(temp_patient_data, new_row)
      
      if(ifPRIASBiopsy(temp_patient_data, visit_schedule[i], previous_biopsy_time)){
        proposed_biopsy_times = c(proposed_biopsy_times, visit_schedule[i])
        previous_biopsy_time = visit_schedule[i]
      }
    }
    
    proposed_biopsy_times_list[[sim]] = proposed_biopsy_times
  }
  
  mode_schedule_string = names(table(unlist(lapply(proposed_biopsy_times_list, paste, collapse = "-"))))[1]
  
  final_proposed_biopsy_times = as.numeric(strsplit(mode_schedule_string, split = '-')[[1]])
  return(final_proposed_biopsy_times)
}

ifFixedRiskBasedBiopsy =  function(object, patient_data, cur_visit_time, 
                                   last_biopsy_time, min_biopsy_gap=1, threshold, M=500){
  
  if(cur_visit_time - last_biopsy_time>=min_biopsy_gap){
    visit_cum_risk = rowMeans(1-getExpectedFutureOutcomes(object, patient_data, last_biopsy_time, 
                                                          survival_predict_times = cur_visit_time,
                                                          psaDist = "Tdist", M = M)$predicted_surv_prob)
    return(visit_cum_risk >= threshold)
  }else{
    return(FALSE)
  }
}

getFixedRiskBasedScheduleNoGrid = function(object, patient_data, 
                                           cur_visit_time, 
                                           last_biopsy_time,
                                           risk_threshold = 0.05, 
                                           min_biopsy_gap = 1,
                                           M=500, horizon=10){
  previous_biopsy_time = last_biopsy_time
  visit_time = cur_visit_time
  proposed_biopsy_times = c()
  while(horizon - previous_biopsy_time >= min_biopsy_gap){
    if(visit_time - previous_biopsy_time >= min_biopsy_gap){
      grid_times = seq(visit_time, horizon, by=0.01)
      visit_cum_risk = rowMeans(1-getExpectedFutureOutcomes(object, patient_data, previous_biopsy_time, 
                                                            survival_predict_times = grid_times,
                                                            psaDist = "Tdist", M = M)$predicted_surv_prob)
      
      nearest_index = which.min(abs(visit_cum_risk - risk_threshold))
      previous_biopsy_time = grid_times[nearest_index]
      
      proposed_biopsy_times = c(proposed_biopsy_times, previous_biopsy_time)
    }
    
    visit_time = previous_biopsy_time + min_biopsy_gap
  }
  
  return(proposed_biopsy_times)
}

#################### CACHE BASED SOLUTION #######################
ifAutomaticRiskBasedBiopsy = function(object, patient_data, cur_visit_time,
                                      last_biopsy_time, min_biopsy_gap = 1, M = M, 
                                      delay_limit = 1, horizon=10, use_exact_delay=T){
  if(cur_visit_time - last_biopsy_time < min_biopsy_gap){
    return(FALSE)
  }
  
  PSA_CHECK_UP_TIME = c(seq(0, 2, 0.25), seq(2.5, horizon, 0.5))
  fixed_grid_visits = c(cur_visit_time, PSA_CHECK_UP_TIME[PSA_CHECK_UP_TIME > cur_visit_time & PSA_CHECK_UP_TIME<=horizon])
  
  post_rand_eff = getRandomEffectDist(object, patient_data, last_biopsy_time, 
                                      psaDist = "Tdist", M = M)
  
  pred_res = getExpectedFutureOutcomes(object, post_rand_eff, patient_data, last_biopsy_time, 
                                       survival_predict_times = fixed_grid_visits,
                                       psaDist = "Tdist", M = M)
  
  fixed_grid_surv = pred_res$predicted_surv_prob
  
  if(use_exact_delay==TRUE){
    fixed_grid_surv = 1 - t(apply(1-fixed_grid_surv, 1, FUN = function(x){
      x / (1-fixed_grid_surv[nrow(fixed_grid_surv),])
    }))
  }
  
  #Fixed risk schedule function
  getRiskScheduleWithCache <- function(surv_threshold){
    proposed_biopsy_times <- c()
    previous_biopsy_time <- last_biopsy_time
    surv_schedule_temp <- fixed_grid_surv
    
    for(i in 1:length(fixed_grid_visits)){
      if(fixed_grid_visits[i] - previous_biopsy_time>=min_biopsy_gap){
        if(mean(surv_schedule_temp[i,], na.rm = T) <= surv_threshold){
          proposed_biopsy_times <- c(proposed_biopsy_times, fixed_grid_visits[i])
          previous_biopsy_time = fixed_grid_visits[i]
          surv_schedule_temp <- apply(fixed_grid_surv, 2, FUN = function(x){x/x[i]})
        }
      }
    }
    
    if(length(proposed_biopsy_times)==0){
      proposed_biopsy_times = horizon
    }else if(horizon - tail(proposed_biopsy_times, 1) < min_biopsy_gap){
      proposed_biopsy_times[length(proposed_biopsy_times)] = horizon
    }else if(tail(proposed_biopsy_times, 1)<horizon){
      proposed_biopsy_times = c(proposed_biopsy_times, horizon)
    }
    
    return(proposed_biopsy_times)
  }
  
  risk_thresholds = seq(0, 1, length.out = 501)
  all_schedules=lapply(1-risk_thresholds, getRiskScheduleWithCache)
  unique_schedules_string = unique(sapply(all_schedules, paste, collapse="-"))
  unique_schedules = lapply(strsplit(unique_schedules_string, split = "-"), as.numeric)
  
  all_biopsy_intervals=lapply(unique_schedules, FUN = function(planned_biopsy_times){
    if(length(planned_biopsy_times)==1){
      biop_intervals = list(c(last_biopsy_time, planned_biopsy_times))
    }else{
      biop_intervals = c(last_biopsy_time,
                         rep(planned_biopsy_times[-length(planned_biopsy_times)],each=2),
                         planned_biopsy_times[length(planned_biopsy_times)])
      biop_intervals = split(biop_intervals, rep(1:(length(biop_intervals)/2), each=2))
    }
    return(biop_intervals)
  })
  unique_biopsy_intervals_string = unique(do.call('c', lapply(all_biopsy_intervals, FUN = function(x){
    sapply(x, paste, collapse="-")
  })))
  
  unique_biopsy_intervals = lapply(strsplit(unique_biopsy_intervals_string, split="-"), as.numeric)
  biopsy_interval_gqpoints = lapply(unique_biopsy_intervals, FUN = getGaussianQuadWeightsPoints)
  biopsy_interval_gqpoints = do.call('c', lapply(biopsy_interval_gqpoints, "[[", "points"))
  
  ###########
  # New surv cache
  ###########
  SURV_CACHE_TIMES = unique(sort(c(last_biopsy_time, fixed_grid_visits, biopsy_interval_gqpoints), decreasing = F))
  pred_res = getExpectedFutureOutcomes(object, post_rand_eff, patient_data, last_biopsy_time, 
                                       survival_predict_times = SURV_CACHE_TIMES[-1],
                                       psaDist = "Tdist", M = M)
  SURV_CACHE_FULL = rbind(rep(1, M), pred_res$predicted_surv_prob)
  
  if(use_exact_delay==TRUE){
    SURV_CACHE_FULL = 1 - t(apply(1-SURV_CACHE_FULL, 1, FUN = function(x){
      x / (1-SURV_CACHE_FULL[nrow(SURV_CACHE_FULL),])
    })) 
  }
  
  #Now we calculate consequences for each unique schedule
  res = vector("list", length=length(unique_schedules))
  names(res) = c(paste("Risk:", unique_schedules))
  for(i in 1:length(unique_schedules)){
    res[[i]] = getConsequencesWithCache(SURV_CACHE_TIMES, SURV_CACHE_FULL, horizon,
                                        unique_schedules[[i]], last_biopsy_time)
  }
  
  expected_delays = sapply(res, "[[", "expected_delay")
  exp_total_biopsies = sapply(res, "[[", "expected_num_biopsies")
  dist = sqrt(expected_delays^2 + (exp_total_biopsies - 1)^2)
  
  # if(length(res)==1){
  #   sd_expected_delay = 1
  #   sd_exp_total_biopsy = 1
  # }else{
  #   sd_expected_delay = sd(expected_delays)
  #   sd_exp_total_biopsy = sd(exp_total_biopsies)
  # }
  
  #scaled_expected_delays = expected_delays/sd_expected_delay
  #scaled_exp_total_biopsies = exp_total_biopsies/sd_exp_total_biopsy
  #scaled_dist = sqrt(scaled_expected_delays^2 + (scaled_exp_total_biopsies - 1/sd_exp_total_biopsy)^2)
  
  # while(sum(scaled_expected_delays<=(delay_limit/sd_expected_delay)) == 0){
  #   delay_limit = delay_limit + 0.1
  # }
  
  #filter = scaled_expected_delays<=(delay_limit/sd_expected_delay)
  #scaled_dist = scaled_dist[filter]
  #res = res[filter]
  # optimal_schedule_index = which.min(scaled_dist)[1]
  
  while(sum(expected_delays<=delay_limit) == 0){
    delay_limit = delay_limit + 0.01
  }
  
  filter = expected_delays<=delay_limit
  dist = dist[filter]
  res = res[filter]
  optimal_schedule_index = which.min(dist)[1]
  
  optimal_schedule = res[[optimal_schedule_index]]$planned_biopsy_times
  decision = optimal_schedule[1] <= cur_visit_time
  
  return(decision)
}

getConsequencesWithCache = function(cache_times, cache_surv_full, horizon,
                                    planned_biopsy_times, latest_survival_time){
  
  #Now we create intervals in which to integrate the conditional surv prob
  if(length(planned_biopsy_times)==1){
    biop_intervals = list(c(latest_survival_time, planned_biopsy_times))
  }else{
    biop_intervals = c(latest_survival_time,
                       rep(planned_biopsy_times[-length(planned_biopsy_times)],each=2),
                       planned_biopsy_times[length(planned_biopsy_times)])
    biop_intervals = split(biop_intervals, rep(1:(length(biop_intervals)/2), each=2))
  }
  
  res = vector("list", length(biop_intervals))
  exp_num_biopsies = rep(0, M)
  for(j in 1:length(biop_intervals)){
    wt_points=getGaussianQuadWeightsPoints(biop_intervals[[j]])
    lower_limit = biop_intervals[[j]][1]
    upper_limit = biop_intervals[[j]][2]
    
    lower_limit_nearest_index = which.min(abs(lower_limit-cache_times))
    upper_limit_nearest_index = which.min(abs(upper_limit-cache_times))
    res[[j]]$cum_risk_interval = cache_surv_full[lower_limit_nearest_index,] - 
      cache_surv_full[upper_limit_nearest_index,]
    
    cond_expected_fail_time = sapply(1:length(wt_points$points), function(i){
      cum_surv_at_points = cache_surv_full[which.min(abs(wt_points$points[i]-cache_times)),] - cache_surv_full[upper_limit_nearest_index,]
      scaled_cum_surv_at_points = cum_surv_at_points/res[[j]]$cum_risk_interval
      
      return(wt_points$weights[i] * scaled_cum_surv_at_points)
    })
    
    res[[j]]$delay = upper_limit - (lower_limit + apply(cond_expected_fail_time, 1, sum))
    exp_num_biopsies = exp_num_biopsies + j * res[[j]]$cum_risk_interval
  }
  exp_num_biopsies = exp_num_biopsies + j * cache_surv_full[upper_limit_nearest_index,]
  exp_num_biopsies = mean(exp_num_biopsies, na.rm = T)
  
  expected_delay = mean(apply(sapply(res, FUN = function(x){
    x$delay * x$cum_risk_interval
  }),1, sum, na.rm=T),na.rm=T)
  
  #proposed biopsy times always has a final biopsy at year 10
  return(list(expected_delay = expected_delay,
              expected_num_biopsies = exp_num_biopsies,
              planned_biopsy_times = planned_biopsy_times))
}

ifAutomaticFullTreeBasedBiopsy = function(object, patient_data, cur_visit_time,
                                          last_biopsy_time, min_biopsy_gap = 1, M = M,
                                          delay_limit = 1, horizon=10, use_exact_delay=T){
  if(cur_visit_time - last_biopsy_time < min_biopsy_gap){
    return(FALSE)
  }
  
  PSA_CHECK_UP_TIME = c(seq(0, 2, 0.25), seq(2.5, horizon, 0.5))
  fixed_grid_visits = c(cur_visit_time, PSA_CHECK_UP_TIME[PSA_CHECK_UP_TIME > cur_visit_time & PSA_CHECK_UP_TIME<=horizon])
  
  #Function to create all possible biopsy schedules
  allBiopsySchedules = function(cur_visit_time, min_biopsy_gap, horizon){
    times = list()
    count = 0
    
    PSA_CHECK_UP_TIME = c(seq(0, 2, 0.25), seq(2.5, horizon, 0.5))
    
    recursive = function(last_biopsy_time, current_visit_time, cur_path){
      if(current_visit_time < cur_visit_time){
        count <<- count + 1
        times[[count]] <<- as.numeric(strsplit(cur_path, split = "-")[[1]])
      }else{
        if(last_biopsy_time - current_visit_time>=min_biopsy_gap){
          #case 1 do a biopsy
          recursive(current_visit_time,
                    current_visit_time - min_biopsy_gap,
                    paste0(current_visit_time,"-",cur_path))
          
          #case 2 don't do a biopsy
          recursive(last_biopsy_time,
                    tail(PSA_CHECK_UP_TIME[PSA_CHECK_UP_TIME<current_visit_time],1),
                    cur_path)
        }
      }
    }
    recursive(horizon, horizon-min_biopsy_gap, as.character(horizon))
    
    #I will sort them by first biopsy time
    first_biopsy_time = sapply(times, min)
    total_biopsies = sapply(times, length)
    times = times[order(first_biopsy_time, total_biopsies, decreasing = F)]
    return(times)
  }
  
  all_schedules = allBiopsySchedules(cur_visit_time = cur_visit_time,
                                    min_biopsy_gap = min_biopsy_gap, horizon = horizon)
  
  unique_schedules_string = unique(sapply(all_schedules, paste, collapse="-"))
  unique_schedules = lapply(strsplit(unique_schedules_string, split = "-"), as.numeric)
  
  all_biopsy_intervals=lapply(unique_schedules, FUN = function(planned_biopsy_times){
    if(length(planned_biopsy_times)==1){
      biop_intervals = list(c(last_biopsy_time, planned_biopsy_times))
    }else{
      biop_intervals = c(last_biopsy_time,
                         rep(planned_biopsy_times[-length(planned_biopsy_times)],each=2),
                         planned_biopsy_times[length(planned_biopsy_times)])
      biop_intervals = split(biop_intervals, rep(1:(length(biop_intervals)/2), each=2))
    }
    return(biop_intervals)
  })
  unique_biopsy_intervals_string = unique(do.call('c', lapply(all_biopsy_intervals, FUN = function(x){
    sapply(x, paste, collapse="-")
  })))
  
  unique_biopsy_intervals = lapply(strsplit(unique_biopsy_intervals_string, split="-"), as.numeric)
  biopsy_interval_gqpoints = lapply(unique_biopsy_intervals, FUN = getGaussianQuadWeightsPoints)
  biopsy_interval_gqpoints = do.call('c', lapply(biopsy_interval_gqpoints, "[[", "points"))
  
  ###########
  # New surv cache
  ###########
  post_rand_eff = getRandomEffectDist(object, patient_data, last_biopsy_time, 
                                      psaDist = "Tdist", M = M)
  
  SURV_CACHE_TIMES = unique(sort(c(last_biopsy_time, fixed_grid_visits, biopsy_interval_gqpoints), decreasing = F))
  pred_res = getExpectedFutureOutcomes(object, post_rand_eff, patient_data, last_biopsy_time, 
                                       survival_predict_times = SURV_CACHE_TIMES[-1],
                                       psaDist = "Tdist", M = M)
  SURV_CACHE_FULL = rbind(rep(1, M), pred_res$predicted_surv_prob)
  
  if(use_exact_delay==TRUE){
    SURV_CACHE_FULL = 1 - t(apply(1-SURV_CACHE_FULL, 1, FUN = function(x){
      x / (1-SURV_CACHE_FULL[nrow(SURV_CACHE_FULL),])
    })) 
  }
  
  #Now we calculate consequences for each unique schedule
  res = vector("list", length=length(unique_schedules))
  names(res) = c(paste("Risk:", unique_schedules))
  for(i in 1:length(unique_schedules)){
    res[[i]] = getConsequencesWithCache(SURV_CACHE_TIMES, SURV_CACHE_FULL, horizon,
                                        unique_schedules[[i]], last_biopsy_time)
  }
  
  expected_delays = sapply(res, "[[", "expected_delay")
  exp_total_biopsies = sapply(res, "[[", "expected_num_biopsies")
  dist = sqrt(expected_delays^2 + (exp_total_biopsies - 1)^2)
  
  # if(length(res)==1){
  #   sd_expected_delay = 1
  #   sd_exp_total_biopsy = 1
  # }else{
  #   sd_expected_delay = sd(expected_delays)
  #   sd_exp_total_biopsy = sd(exp_total_biopsies)
  # }
  
  #scaled_expected_delays = expected_delays/sd_expected_delay
  #scaled_exp_total_biopsies = exp_total_biopsies/sd_exp_total_biopsy
  #scaled_dist = sqrt(scaled_expected_delays^2 + (scaled_exp_total_biopsies - 1/sd_exp_total_biopsy)^2)
  
  # while(sum(scaled_expected_delays<=(delay_limit/sd_expected_delay)) == 0){
  #   delay_limit = delay_limit + 0.1
  # }
  
  #filter = scaled_expected_delays<=(delay_limit/sd_expected_delay)
  #scaled_dist = scaled_dist[filter]
  #res = res[filter]
  # optimal_schedule_index = which.min(scaled_dist)[1]
  
  while(sum(expected_delays<=delay_limit) == 0){
    delay_limit = delay_limit + 0.01
  }
  
  filter = expected_delays<=delay_limit
  dist = dist[filter]
  res = res[filter]
  optimal_schedule_index = which.min(dist)[1]
  
  optimal_schedule = res[[optimal_schedule_index]]$planned_biopsy_times
  decision = optimal_schedule[1] <= cur_visit_time
  
  return(decision)
}