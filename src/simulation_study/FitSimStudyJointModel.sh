#!/bin/bash

#In our paper we did 500 simulations total
total_simulation_rounds=500

#Per simulation we fitted one joint model using 2 cores for each Rscript
max_cores_per_simulation=2

#We cannot run all 500 simulations together, we run 5 of them at a time.
max_simulations_in_parallel=5
total_batches=$((total_simulation_rounds/max_simulations_in_parallel))

#we run simulations with different seeds, starting from 2001
#why 2001? choice of seeds in our paper are motivated by the current year and century
start_seed=2001
end_seed=$((start_seed+total_simulation_rounds-1))

for((batch=0;batch<total_batches;batch++))
  do
    batch_start_seed=$((start_seed + batch*max_simulations_in_parallel))
    batch_end_seed=$((batch_start_seed + max_simulations_in_parallel-1))
    for((seed=batch_start_seed;seed<=batch_end_seed;seed++))
    do
      (
        Rscript FitSimStudyJointModel.R $seed $max_cores_per_simulation > "log_model/simmodel_${seed}.log" 2>&1
      ) &
    done    
    wait
done