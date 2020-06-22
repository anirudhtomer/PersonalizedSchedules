#!/bin/bash
total_test_patients=250
batch_size=5
max_batch_number=$((total_test_patients/batch_size))

seed_start=2001
seed_end=2500

for ((seed=seed_start; seed<=seed_end; seed++))
do
	for ((batch_number=1; batch_number<=max_batch_number; batch_number++))
	do
		for ((pat_counter=1; pat_counter<=batch_size; pat_counter++))
		do
		(
		  testId=$(( 750 + pat_counter + batch_size * (batch_number-1) ))
		  Rscript RunBiopsySchedules.R $seed $testId > "log_schedules/pat_${seed}_${testId}.log" 2>&1
		) &
		done
		wait
	done
done

echo "Finished all simulations. Combining results now."
Rscript CombinedResults.R $seed_start $seed_end