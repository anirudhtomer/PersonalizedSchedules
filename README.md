# Personalized Schedules for Burdensome Surveillance Tests
Code and (synthetic) data accompanying manuscript titled 'Personalized Schedules for Shared Decision Making of Burdensome Surveillance Tests'

# Author Contributions Checklist Form
## Data

### Abstract
In this work, we developed personalized test schedules for invasive diagnostic tests in the surveillance of chronic non-communicable diseases. We demonstrate our methodology for the prostate cancer surveillance scenario, using a model fitted to the world's largest prostate cancer active surveillance dataset abbreviated as PRIAS, dated April 2019. More than 100 medical centers from 17 countries contribute to PRIAS using a common study protocol, and web-based tool, both available at www.prias-project.org. The PRIAS dataset consists of patient age at inclusion, and repeated measurements data on biopsy Gleason grade groups (invasive diagnostic test), prostate-specific antigen or PSA (continuous: biomarker ng/mL), digital rectal examination or DRE (binary: tumor palpable or not), and magnetic resonance imaging or MRI (PI-RADS score).

### Availability
The PRIAS dataset is not openly accessible, but access to the dataset can be requested based on a study proposal approved by the PRIAS steering committee. The website of the PRIAS program is www.prias-project.org. We do provide a synthetic dataset instead of the original dataset. Readers can use this dataset along with the provided code to understand and recreate the workflow required to make personalized schedules in their datasets.

The **main contribution of our work** is not the model fitted to the PRIAS dataset, but rather the methodology proposed to schedule invasive diagnostic tests. The methodology is generic for use in other surveillance scenarios (e.g., scheduling endoscopies in Barett's esophagus). Hence, the lack of direct public access to the prostate cancer dataset PRIAS, which we only used for demonstration purposes, is not a major limitation. The simulation study we conducted utilizes the parameter estimates of the model fitted to the PRIAS dataset and can be reproduced without access to the PRIAS dataset itself.

### Description
Upon inclusion in PRIAS, patients undergo repeat biopsies at year one, four, seven, and ten of follow-up. Additional yearly biopsies are scheduled when a patient's PSA doubling time is between zero and ten years. Biopsy results are measured as Gleason grade groups (1 to 5). An inclusion criterion in PRIAS is that patient's Gleason grade group should be equal to 1. We selected all 7813 patients who had Gleason grade group 1 at inclusion. Our primary event of interest is an increase in this Gleason grade group observed upon repeat biopsy, called progression (1134 patients). Progression is a trigger for treatment advice in PRIAS. In PRIAS 2250 patients were provided treatment based on their PSA, the number of biopsy cores with cancer, or anxiety/other reasons. However, our reasons for focusing solely on progression are that progression is strongly associated with cancer-related outcomes, and other treatment triggers vary between cohorts.

The PSA is measured every three months for the first two years and every six months after that. The DRE is conducted every six months. We use PSA after using log<sub>2</sub>(PSA + 1) transformation (details in the manuscript). Due to many sparse categories of the ordinal DRE data, we use DRE as a binary variable, namely as DRE=T1c versus DRE > T1c. A DRE measurement equal to stage T1c indicates a clinically inapparent tumor that is not palpable or visible by imaging. In contrast, tumors with  DRE > T1c are palpable. We excluded MRI data because it is available very sparsely.

## Code for Reproducing Results and Creating Personalized Test Schedules
We provide R code for the following purposes:
1. To generate synthetic PRIAS like datasets.
2. For cleaning the original PRIAS dataset for readers who obtain a copy of the raw dataset from the PRIAS consortium.
3. To fit a joint model to the synthetic or real PRIAS dataset using the R package **[JMBayes](https://cran.r-project.org/web/packages/JMbayes/index.html)**. 
4. For reproducing the simulation study results comparing different biopsy schedules in PRIAS. 
5. A generic API that can be used to schedule invasive tests in a personalized manner for any surveillance scenario. This API is compatible with R joint model objects fitted using the JMbayes package. 
6. Code to reproduce figures and tables.

### Must Do's Before Running The Code
This repository has been divided into two main folders, namely [data](data) and [src](src), for Rdata files and source code, respectively. Readers must clone the complete repository and then set the downloaded folder's location as the current working directory in R. For example,
```
$ setwd('/home/username_here/PersonalizedSchedules')
```
**This is also the first command in each of the source code files. Hence, readers must update the working directory in each soure code files for it to work.**

In total the provided source code depends on four R packages:
- [JMBayes](https://cran.r-project.org/web/packages/JMbayes/index.html)
- [survival](https://cran.r-project.org/web/packages/survival/index.html)
- [MASS](https://cran.r-project.org/web/packages/MASS/index.html)
- splines: part of the R core now.

The code has been tested to work with 64 bit R-3.6.1 installed on a laptop with 4th generation Intel Core-i7 processor, and 8GB RAM.

### Files That Should Not Be Changed or Removed
#### Joint model fitted to the original PRIAS dataset
The file [data/fitted_model_original_prias_dataset.Rdata](data/fitted_model_original_prias_dataset.Rdata) contains the MCMC samples for each of the parameter (excluding random effects) in the joint model fitted to the PRIAS dataset. For example, readers can check the mean and 95% confidence interval of the posterior predictive distribution of the parameters using the R command:
```
> load('data/fitted_model_original_prias_dataset.Rdata')
> lapply(fitted_model_original_prias_dataset$mcmc, FUN = mean)
> lapply(fitted_model_original_prias_dataset$mcmc, FUN = quantile, probs=c(0.025, 0.975))
```
This fitted joint model is used to generate synthetic PRIAS like datasets, as well as to generate simulation study datasets. Hence, we advise not to replace or delete this file.

#### Results of the simulation study we conducted for 500 x 250 test patients
The results of the simulation study we conducted are saved in the file [simstudy_results.Rdata](data/simulation_study/simstudy_results.Rdata). The dataframe object contained in this file is called *simstudy_results* and it consists of 7 columns. These columns are:
- *seed* is the seed of with which the whole test dataset of 250 patients was generated.
- *P_ID* is patient ID (between 751 and 1000) in the dataset of 250 patients.
- *progression_time* is the true time of progression (simulated) of the test patient.
- *progressed* is a binary indicator showing if the patient has progression time less than the ten year follow-up period for which we conduct simulation study.
- *schedule* is the name of a schedule and it can take following values: Annual (yearly biopsy), PRIAS (PRIAS schedule), FixedRisk_5% (risk based biopsy with fixed 5% threshold), FixedRisk_10% (risk based biopsy with fixed 10% threshold), AutoRisk_Delay_0.75 (risk based biopsy with automatically chosen threshold based on our utility function, although with the constraint that average delay is less than 0.75 years), and AutoRisk_NoDelayLimit (risk based biopsy with automatically chosen threshold based on our utility function).
- *nb* is the number of biopsies conducted.
- *delay* is the time delay in detecting progression, and isdefined only for patients who progressed.

#### Utility files
Our source code relies on the following utility files:
- [VariousSchedules.R](src/utility_files/VariousSchedules.R) contains code for the following schedules:
    - Fixed: annual, and PRIAS biopsy schedules. 
    - Risk based: risk based schedules with fixed threshold, and risk based schedules based on a dynamically selected risk threshold using the utility function proposed in our manuscript.

- [SyntheticDataSettingsAndFunctions.R](src/utility_files/SyntheticDataSettingsAndFunctions.R) contains code for generating sythentic data, and fitting joiint model to the synthetic datasets.
- [prediction_psa_dre.R](src/utility_files/prediction_psa_dre.R) and [prediction_psa_dre_randEff_reuse.R](src/utility_files/prediction_psa_dre_randEff_reuse.R) contain code for predicting cumulative-risk of disease progression and prostate-specific antigen given a patient's observed longitudinal data.

<hr>

### 1. Generating synthetic PRIAS like synthetic datasets
Although the PRIAS dataset is not publicly available, PRIAS like synthetic dataset can be generated using the R code given in the file [GenerateSyntheticData.R](src/data_analysis/GenerateSyntheticData.R). Within this file the loss to follow-up process can be controlled using the variable *loss_to_follow_up_max_time* on line 39. While we assume a uniform distribution for loss to follow-up, readers can change this distribution on line 42. Loss to follow-up before year one of follow-up is very rare in the prostate cancer surveillance scenario. On the contrary it is mandtory to have a confirmatory biopsy at year one. Hence loss to follow-up times should be only after year one. To generate a synthetic dataset of 1000 patients the following command can be used. Running this command will take 3-4 minutes approximately.
```
$ Rscript src/data_analysis/GenerateSyntheticData.R 1000
```
Upon completion, a synthetic patient R dataframe will be saved with the name *longitudinal_data* in the so called long-format in the file [data/synthetic_data.Rdata](data/synthetic_data.Rdata). By default, we have uploaded a synthetic dataset of 1000 patients. The data of first patient in this dataset is shown below. The data consists of the following 7 variables:
- *P_ID* is patient ID.
- *age* is patient age at inclusion in surveillance.
- *year_visit* is the time of the follow-up visit at which PSA or DRE is measured.
- *log2psaplus1* are log<sub>2</sub>(PSA + 1) transformed PSA values.
- *palpable_dre* is a binary variable indicating whether tumor at a visit was palpable or not.
- *latest_survival_time* and *earliest_failure_time* indicate the interval-censored time of disease progression of the patient. This interval censored time is generated using the PRIAS biopsy schedule, implemented in the function *ifPRIASBiopsy* in the file [VariousSchedules.R](src/utility_files/VariousSchedules.R).

Example dataset for one synthetic patient is shown below.

| P_ID | age   | year_visit | log2psaplus1 | palpable_dre | latest_survival_time | earliest_failure_time |
|------|-------|------------|--------------|--------------|----------------------|-----------------------|
| 1    | 71.92 | 0          | 3.13         | 0            | 5                    | Inf                   |
| 1    | 71.92 | 0.25       | 2.7          | NA           | 5                    | Inf                   |
| 1    | 71.92 | 0.5        | 3.09         | 0            | 5                    | Inf                   |
| 1    | 71.92 | 0.75       | 3.38         | NA           | 5                    | Inf                   |
| 1    | 71.92 | 1          | 3.32         | 0            | 5                    | Inf                   |
| 1    | 71.92 | 1.25       | 3.22         | NA           | 5                    | Inf                   |
| 1    | 71.92 | 1.5        | 3.34         | 0            | 5                    | Inf                   |
| 1    | 71.92 | 1.75       | 3.39         | NA           | 5                    | Inf                   |
| 1    | 71.92 | 2          | 3.34         | 0            | 5                    | Inf                   |
| 1    | 71.92 | 2.5        | 3.39         | 0            | 5                    | Inf                   |
| 1    | 71.92 | 3          | 3.87         | 0            | 5                    | Inf                   |
| 1    | 71.92 | 3.5        | 3.95         | 0            | 5                    | Inf                   |
| 1    | 71.92 | 4          | 3.83         | 0            | 5                    | Inf                   |
| 1    | 71.92 | 4.5        | 3.95         | 0            | 5                    | Inf                   |
| 1    | 71.92 | 5          | 3.87         | 0            | 5                    | Inf                   |

### 2. Data Cleaning For the Real PRIAS Dataset
If readers request the PRIAS dataset (April 2019 version), it may be delivered in SPSS's SAV format, with one data row per patient. To clean this dataset and to convert it in the long format required to fit a joint model to it please use the file [PRIAS_SPSSDataCleaningAndFormatting.R](src/data_analysis/PRIAS_SPSSDataCleaningAndFormatting.R). Suppose the raw PRIAS dataset is saved with the name `raw_prias_data.sav' in the [data](data) folder. Then the following command will clean this dataset:
```
$ Rscript src/data_analysis/PRIAS_SPSSDataCleaningAndFormatting.R data/raw_prias_data.sav
```
The cleaned dataset is named *longitudinal_data*, and will be saved in the long format at the location [data/clean_longitudinal_data.Rdata](data/clean_longitudinal_data.Rdata). Currently there is no dataset at this location, and so the link does not work.

### 3. Analysis of Synthetic and/or real PRIAS Dataset
The analysis of synthetic or real PRIAS dataset using a joint model for time-to-event and longitudinal outcomes is done in [AnalysisSyntheticData.R](src/data_analysis/AnalysisSyntheticData.R). The file expects that the patient dataframe is named *longitudinal_data* and is in the long format. It is important to note that it takes a laptop computer with the specification mentioned above between 8-10 hours to fit a model to the real PRIAS dataset using the R package JMbayes. For a synthetic dataset with 1000 patients it takes between 2 to 3 hours. To fit a joint model to the uploaded [synthetic dataset](data/synthetic_data.Rdata) of 1000 patients the following command can be used:
```
$ Rscript src/data_analysis/AnalysisSyntheticData.R data/synthetic_data.Rdata
```
The fitted joint model object (see [JMbayes]((https://cran.r-project.org/web/packages/JMbayes/index.html))) is given the name *fitted_jointmodel* and saved at the location [data/fitted_jointmodel.Rdata](data/fitted_jointmodel.Rdata).

### 4. Reproducing Simulation Study Results
The simulation study we conducted in our work required running it on Linux server computers. We conducted 500 simulations, each with total 1000 patients: 750 training patients and 250 test patients. The simuation study consists of two sequential steps. First, to generate simulated data and fitting a joint model to the training patients. Second, to use the fitted joint model to conduct conduct a randomized clinical trial of different fixed and personalized biopsy schedules for test patients. 

#### Step 4.1: Generating 500 simulated datasets and fitting joint models to training patients
For this purpose readers can utilize the script file [FitSimStudyJointModel.sh](src/simulation_study/FitSimStudyJointModel.sh). Example:
```
$ chmod u+x FitSimStudyJointModel.sh
$ ./FitSimStudyJointModel.sh
```
The shell script code uses seeds between 2001 and 2500 for the 500 simulated datasets. The seed range can be changed using the variables *start_seed* and *total_simulation_rounds*, currently set to 2001 and 500 respectively. Internally the shell script invokes the R code in the file [FitSimStudyJointModel.R](src/simulation_study/FitSimStudyJointModel.R). This file generates a dataset of 1000 patients with a particular seed, divides the dataset into training (750 patients) and test (250) parts, fits a joint model to the training patients using R package JMbayes, and finally saves the simulation data with the object name *jointModelData* in the [data/simulation_study/fitted_models/](data/simulation_study/fitted_models/) folder. Each saved file is uniquely identifiable with the seed with the saved data was generated. For example, simulation data generated using a seed 2001 is saved as 'jointModelData_seed_2001.Rdata'. During the simulation all logs are saved in the folder [src/simulation_study/log_model](src/simulation_study/log_model), with one log file, e.g., 'simmodel_2001.log', per simulated dataset.

The **memory requirement** per simulation is around 1.2GB. The joint model objects are also quite large in size, and while saving take around 300MB of space on the hard disk. To conserve hard disk space, before saving joint model objects we trim them. However, after trimming, the objects do not work with the JMbayes package anymore. The trimming code is between line 39 and line 49 in the file [FitSimStudyJointModel.R](src/simulation_study/FitSimStudyJointModel.R). Commenting out these lines does not affect results or change workflow, but they do increase memry requirements for the step 2 of the simulation study where we use fitted models to run biopsy schedules. Due to the memory requirements most servers cannot run all 500 simulations together. Hence the script [FitSimStudyJointModel.sh](src/simulation_study/FitSimStudyJointModel.sh) allows running simulations in batches. The batch size (default 5) can be controlled via the variable *max_simulations_in_parallel*. Please note that each simulation requires minimum of 2 CPU cores (on Linux 2 R processes). That is, on a Linux server, a batch of 5 simulations together will lead to 10 R processes in total. It takes 75 minutes to generate simulation data and to fit a joint model to it with each CPU core having a maximum clock speed of 2GHz, and total server RAM of 65GB.

#### Step 4.2: Using the fitted joint models to conduct a randomized clinical trial of biopsy schedules for test patients
In this step we conduct a randomized clinical trial of biopsy schedules for each of the 500 x 250 test patients generated in previous step. The different schedules we employ are, namely, annual schedule, PRIAS schedule, risk based schedule with 5% and 10% thresholds, and risk based schedules based on the utility function given in our manuscript. In this regard, the file [RunBiopsySchedules.R](src/simulation_study/RunBiopsySchedules.R) runs all of these schedules for a single patient. Each of the 250 test patient patients has an ID between 751 and 1000, and is associated with a particular seed between 2001 and 2500 using which the test patients were generated. That is, the unique identifier of a patient is the seed and the patient ID. Running a simulation study, for a patient with ID 999 from the dataset with seed 2001 can be done with the command:
```
$ Rscript RunBiopsySchedules.R 2001 999
```
For the selected patient, this R code saves timings of biopsies in a dataframe with object name *biopsyDf*. The dataframe is saved in the folder [data/simulation_study/schedule_results/](data/simulation_study/schedule_results/), with filename uniquely identifiable, e.g., 'patient_seed_2001_P_ID_999.Rdata'. During the simulation all logs are saved in the folder [src/simulation_study/log_schedules](src/simulation_study/log_schedules), with one log file, e.g., 'pat_2001_999.log', per patient.

Schedules can be run for all 500 x 250 patients together using the shell script [RunBiopsySchedules.sh](src/simulation_study/RunBiopsySchedules.sh). For example:
```
$ chmod u+x RunBiopsySchedules.sh
$ ./RunBiopsySchedules.sh
```
By default the shell script runs schedules for a batch of 5 patients at a time, with one R process for each patient's schedules. This batch size can be changed in the variable *batch_size* in the file [RunBiopsySchedules.sh](src/simulation_study/RunBiopsySchedules.sh). It takes 10 minutes to conduct all schedules for a single patient on a single CPU core having a maximum clock speed of 2GHz, and total server RAM of 65GB. That is, nearly 868 days for a complete simulation study of 500 x 250 patients. We reduced this time to one month by running simulations on a supercomputer using the provided shell script.

When all schedules for all patients are finished the shell script automatically combines and summarizes individual patient results saved in [data/simulation_study/schedule_results/](data/simulation_study/schedule_results/) into a single R dataframe, and saves it in a file titled [simstudy_results.Rdata](data/simulation_study/simstudy_results.Rdata). The dataframe object contained in this file is called *simstudy_results* and it consists of 7 columns. These columns are [defined earlier](#results-of-the-simulation-study-we-conducted-for-500-x-250-test-patients).

### 5. API's for Creating Personalized Schedules for Any Surveillance Scenario
We provide two API's for creating personalized schedules for burdensome diagnostic tests in any surveillance scenario and not just prostate cancer. These API's are compatible (and dependent) with joint model objects fitting using the *mvJointModelBayes* function of R package [JMBayes](https://cran.r-project.org/web/packages/JMbayes/index.html). The two API's are, namely:
- *getBurdenAndBenefit* in the file [GetBurdenAndBenefit.R](src/api/GetBurdenAndBenefit.R).
    - Description: This API returns the expected burden and expected benefit of any given planned schedule of tests (fixed or personalized), under the assumption that progression happens at or before the time of the final planned tests. The expected burden is quantitied by the expected number of tests, and expected benefit is quantified by the expected time delay in detecting progression.
    - Arguments:
        - *object*: An object of the class mvJointModelBayes from the [JMBayes](https://cran.r-project.org/web/packages/JMbayes/index.html) package.
        - *newdata*: R dataframe containing data of a new patient for whom we intend to calculate burden and benefit. It should be sorted in ascending order by time of visit.
        - *idVar (deafult='id')*: the name of the column which contains patient's ID in the *newdata*.
        - *last_test_time (default=NULL)*: the time of the last invasive test at which progression was not detected. When set to NULL, last_test_time is set to the time of the last visit of the patient.
        - *planned_test_schedule*: the test schedule for which we intend to calculate expected time delay in detecting progression (shorter is beneficial) and expected number of invasive tests (burden).
        - *seed (default=1L)*: numeric scalar, the random seed used to produce the results. see R documentaion for survfitJM in [JMBayes](https://cran.r-project.org/web/packages/JMbayes/index.html) for details.
        - *M (default=400L)*: integer denoting how many Monte Carlo samples to use. see R documentaion for survfitJM in [JMBayes](https://cran.r-project.org/web/packages/JMbayes/index.html) for details.
        - *cache_size (default=1000)*: controls the accuracy and speed of calculations. Internally the size of the cache determines the number of time points between *last_test_time* and *max(planned_test_schedule)* on which we pre-calculate cumulative-risk of disease progression.   Higher cache_size provides more accurate results, but also leads to a much higher computational time.
    - Return Value (list):
        - *expected_detection_delay*: expected time delay in detecting progression if patient progresses at or before *max(planned_test_schedule)*.
        - *expected_num_tests*: expected number of tests that will be conducted until the time *max(planned_test_schedule)*.
        - *planned_test_schedule*: same as the argument that was passed to the function.
        - *euclidean_distance*: Euclidean distance between the points *(1,0)* and *(expected_num_tests, expected_detection_delay)* in a two-dimensional space of number of tests and time delay in detecting progression.
        
- *getPersonalizedSchedule* in the file [GetPersonalizedSchedule.R](src/api/GetPersonalizedSchedule.R).
    - Description: This function creates a personalized schedule of invasive tests using the utility function proposed in our manuscript.
    - Arguments:
        - *object*: An object of the class mvJointModelBayes from the [JMBayes](https://cran.r-project.org/web/packages/JMbayes/index.html) package.
        - *newdata*: R dataframe containing data of a new patient for whom we intend to calculate burden and benefit. It should be sorted in ascending order by time of visit.
        - *idVar (deafult='id')*: the name of the column which contains patient's ID in the *newdata*.
        - *last_test_time (default=NULL)*: the time of the last invasive test at which progression was not detected. When set to NULL, last_test_time is set to the time of the last visit of the patient.
        - *fixed_grid_visits*: the fixed_schedule (*U* in our manuscript) according to which patient will visit for measurement of biomarkers. We decide for invasive tests on each of the *fixed_grid_visits* using a risk based approach.
        - *gap (default=0)*: required minimum time gap between consecutive tests. For example, in prostate cancer it is recommended that tests should be planned with a one year gap between them.
        - *detection_delay_limit (default=Infinity)*: should the Euclidean distance be minimized with a constraint on the expected delay in detecting progression. In other words this is the maximum acceptable limit for expected time delay in detecting progression.
        - *total_tests_limit (default=Infinity)*: should the Euclidean distance be minimized with a constraint on the expected number of tests. In other words this is the maximum acceptable limit for expected number of tests that will be conducted.
        - *seed (default=1L)*: numeric scalar, the random seed used to produce the results. see R documentaion for survfitJM in [JMBayes](https://cran.r-project.org/web/packages/JMbayes/index.html) for details.
        - *M (default=400L)*: integer denoting how many Monte Carlo samples to use. see R documentaion for survfitJM in [JMBayes](https://cran.r-project.org/web/packages/JMbayes/index.html) for details.
        - *cache_size (default=1000)*: controls the accuracy and speed of calculations. Internally the size of the cache determines the number of time points between *last_test_time* and *max(planned_test_schedule)* on which we pre-calculate cumulative-risk of disease progression.   Higher cache_size provides more accurate results, but also leads to a much higher computational time.
    - Return Value (list):
        - *all_schedules*: a list of all possible risk based schedules that we evaluated. Each member of the list contains *expected_detection_delay*, *expected_num_tests*, and *planned_test_schedule* These have been already defined for the previous API titled *getBurdenAndBenefit*.
        - *optimal_schedule*: schedule which optimizes the Euclidean distance between point (1,0) and (*expected_num_tests*, *expected_detection_delay*) among *all_schedules* in a two-dimensional space of number of tests and time delay in detecting progression.
        
### 6. Generating figures and tables presented in our work
Figure 1 to Figure 4 in the main manuscript are only illustrative, and hence we did not provide code for them. 
- Figure 5 can be reproduced using the code [DemoSchedulePlot.R](src/data_analysis/DemoSchedulePlot.R). The plot requires the original PRIAS based fitted joint model. This joint model in its [trimmed form](#joint-model-fitted-to-the-original-prias-dataset) has been uploaded with this repository. However, the trimmed version is not enough to reproduce the figure. We could not upload the untrimmed version because it contains the original PRIAS dataset. An alternative is to fit a joint model to a synthetic dataset using the [AnalysisSyntheticData.R](src/data_analysis/AnalysisSyntheticData.R) file. This file will save a joint model at the location [data/fitted_jointmodel.Rdata](data/fitted_jointmodel.Rdata). The [DemoSchedulePlot.R](src/data_analysis/DemoSchedulePlot.R) requires a joint model at the location [data/fitted_jointmodel.Rdata](data/fitted_jointmodel.Rdata).
- Figure 6 can be reproduced using the code [PlotResults.R](src/simulation_study/PlotResults.R). This code requires on R packages ggplot, and ggpubr. The code also expects simulation study results to be at the location [data/simulation_study/simstudy_results.Rdata](data/simulation_study/simstudy_results.Rdata) in the format [defined earlier](#results-of-the-simulation-study-we-conducted-for-500-x-250-test-patients).
- Supplementary Figure 1 can be reproduced using the file [NPMLEPlot.R](src/data_analysis/NPMLEPlot.R). Because the figure relies on synthetic data and not actual PRIAS dataset, the resulting figure only looks similar to the Figure 1 in the manuscript.
- Supplementary tables pertaining to parameter estimates can be verified using the [code provided above](#joint-model-fitted-to-the-original-prias-dataset).
- Supplementary Table 8 for area under receiver operating characteristic curve (AUC) and mean absolute prediction error (MAPE) estimates can be reproduced using the APIs for AUC and MAPE provided by the [JMBayes](https://cran.r-project.org/web/packages/JMbayes/index.html) package. Since PRIAS data is not uploaded with this manuscript, we advise readers to fit a model to the synthetic dataset and calculate AUC and MAPE using the synthetic dataset and the fitted model.
- Supplmentary Table 9 for simulation results can be verified using the following R code:

```
> load('data/simulation_study/simstudy_results.Rdata')
> # summary of number of biopsies for non-progressing patients
> by(INDICES = simstudy_results$schedule[simstudy_results$progressed==0],
   data= simstudy_results$nb[simstudy_results$progressed==0], summary)
> # summary of number of biopsies for progressing patients
> by(INDICES = simstudy_results$schedule[simstudy_results$progressed==1],
   data= simstudy_results$nb[simstudy_results$progressed==1], summary)
> # summary of time delay in detecting progression for progressing patients
> by(INDICES = simstudy_results$schedule[simstudy_results$progressed==1],
   data= simstudy_results$delay[simstudy_results$progressed==1], summary)
```
