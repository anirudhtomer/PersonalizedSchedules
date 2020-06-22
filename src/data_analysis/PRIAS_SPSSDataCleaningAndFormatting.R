setwd('/home/username_here/PersonalizedSchedules')
if(!require(foreign)){
  install.packages("foreign")
}

rm(list=ls())

args = commandArgs(trailingOnly=TRUE)
prias=read.spss(file = args[1], to.data.frame = T)
YEAR_DIVIDER = 24 * 60 * 60 * 365

prias_mri = prias[,c(1,which(grepl(colnames(prias), pattern = 'Pirad', fixed = F)==T))]
prias_mri$N_A = NA

usefulCols = c(1,2,3,4,6,7,8,14,15,16, 21,22,23,24,25,32,33,34,35,36,38,
                   63:200, 201:292, 385:430, 707:844,
                   1075:1120, 1213:1258, 1949:1994, 1998,1999, 2000, 2001,
                   2002, 2003, 2021, 2022, 2023)

colnames(prias[,usefulCols])

prias = prias[,usefulCols]
prias$Date_birth = (prias$Date_dianosis - prias$Date_birth)/YEAR_DIVIDER
colnames(prias)[2] = "age"
colnames(prias)[11] = "Treatment_Type"
colnames(prias)[12] = "Treatment_Reason"
levels(prias$Treatment_Reason) = trimws(levels(prias$Treatment_Reason))
levels(prias$Treatment_Type) = trimws(levels(prias$Treatment_Type))

#Check how many NA
table(prias$Treatment_Reason, useNA = 'always')
prias$Treatment_Reason[prias$Treatment_Reason %in% c("", "N/A")] = NA
prias$Treatment_Reason=droplevels(prias$Treatment_Reason)

table(prias$Treatment_Type, useNA = 'always')
prias$Treatment_Type[prias$Treatment_Type %in% c("", "N/A")] = NA
prias$Treatment_Type=droplevels(prias$Treatment_Type)

#Do not change the order
prias$Gleason_sum_2 = prias$Gleason1_2 + prias$Gleason2_2
prias$N_A = NA
prias$gleason_rad_prost = prias$Gleason1_Rad_Prost + prias$Gleason2_Rad_Prost

#Remove patients with NA, and age less than 30 and age more than 90 (includes negative ages)
prias = prias[!is.na(prias$age) & prias$age >= 30 & prias$age<=90,]

#16 patients had a repeat Gleason measurement
gl2_filter = !is.na(prias$Gleason_sum_2) & !is.na(prias$Date_dianosis2) & 
  prias$Gleason_sum_2 > 0
sum(gl2_filter)

#lets check these 16 patients Gleason first and second time
prias$Gleason_sum[gl2_filter]
prias$Gleason_sum_2[gl2_filter]
prias$Num_cores2[gl2_filter]
prias$Num_Cores_PC2[gl2_filter]
(prias$Date_dianosis2[gl2_filter] - prias$Date_dianosis[gl2_filter])/YEAR_DIVIDER
prias$P_ID[gl2_filter]

#No patient had a Gleason downgrade upon rebiopsy, but it seems that
#the 4th patient in this list should never have been there in the first place
prias = prias[prias$P_ID != 317,]

#No need to worry about Gleason_sum_2 again

#Lets remove all patients who have a Gleason > 6 at first visit
prias = prias[!(prias$Gleason_sum %in% c(7,8,9,10)),]

#The column 584 is a dummy column to handle Gleason2 which is 
#gleason taken again at the time of diagnosis
prias_long=reshape(prias, direction='long', idvar='P_ID', timevar = "visit_number",
                   varying=list(c(3, 584, 68:113), c(4, 584, 22:67), 
                                c(3, 18, 252:297), c(5, 584, 114:159), c(5, 584, 528:573),
                                c(8, 19, 298:343), c(9, 20, 344:389), c(10, 583, 390:435), 
                                c(6, 16, 160:205), c(7, 17, 206:251), 
                                c(584, 584, 436:481), c(21, 584, 482:527)),
                   v.names=c('dom_psa', 'psa', 
                             'dom_dre_gleason', 'dre', 'dre_recode',
                             'gleason_1', 'gleason_2', 'gleason_sum', 
                             'ncores', 'ncores_pc', 
                             'dummy', 'active'))

prias_long = prias_long[order(prias_long$P_ID), ]

prias_mri_long = reshape(prias_mri, direction='long', idvar='P_ID', timevar = "visit_number",
                         varying=list(c(2,143,3:48), c(49, 143, 50:95), 
                                      c(96, 143, 97:142)), 
                         v.names=c('pirads1', 'pirads2','pirads3'))
prias_mri_long = prias_mri_long[order(prias_mri_long$P_ID), ]
prias_mri_long = prias_mri_long[prias_mri_long$P_ID %in% prias$P_ID,]

#prias_long = cbind(prias_long, prias_mri_long[, c('pirads1', 'pirads2','pirads3')])

#First lets check Gleason table
table(prias_long$gleason_sum, useNA = "always")

#Didn't find anything negative. Now lets check ncores and ncores_pc
table(prias_long$ncores, useNA = "always")
table(prias_long$ncores_pc, useNA = "always")

# No one entered ncores_pc=0 as NA. thank god
sum(!is.na(prias_long$ncores) & is.na(prias_long$ncores_pc))

#Some cores are NA, 99, 999 or 0, and some cores_pc are 99 or NA. 
#Lets set those Gleason and cores as NA
temp_filter = prias_long$ncores %in% c(NA, 99, 999, 0) | prias_long$ncores_pc %in% c(NA, 99)
prias_long$gleason_1[temp_filter] = NA
prias_long$gleason_2[temp_filter] = NA
prias_long$gleason_sum[temp_filter] = NA
prias_long$ncores[temp_filter] = NA
prias_long$ncores_pc[temp_filter] = NA

#Sometimes ncores with pc are more than total number of cores taken.
#Lets set those Gleason and cores as NA
summary(prias_long$ncores - prias_long$ncores_pc)

temp_filter = is.na(prias_long$ncores - prias_long$ncores_pc) | (prias_long$ncores - prias_long$ncores_pc) < 0
prias_long$gleason_1[temp_filter] = NA
prias_long$gleason_2[temp_filter] = NA
prias_long$gleason_sum[temp_filter] = NA
prias_long$ncores[temp_filter] = NA
prias_long$ncores_pc[temp_filter] = NA

#Lets check first row data, we have 48 rows per patient
firstRowLongFilter = (prias_long$visit_number==1)

#Set Gleason = 6 for all patients on row 1
prias_long$gleason_sum[firstRowLongFilter] = 6

# add a date of diagnosis in long
# date of first psa and first dre gleason is same
summary(prias_long$dom_dre_gleason[firstRowLongFilter] -prias_long$dom_psa[firstRowLongFilter])
prias_long$dom_diagnosis = rep(prias_long$dom_dre_gleason[firstRowLongFilter], each=48)

prias_long$year_dre_gleason = (prias_long$dom_dre_gleason - prias_long$dom_diagnosis)/YEAR_DIVIDER
prias_long$year_psa = (prias_long$dom_psa - prias_long$dom_diagnosis)/YEAR_DIVIDER
prias_long$year_surgery = (prias_long$Surgerydate - prias_long$dom_diagnosis)/YEAR_DIVIDER
prias_long$year_dead = (prias_long$Dead_Date - prias_long$dom_diagnosis)/YEAR_DIVIDER
prias_long$year_discontinued = (prias_long$Date_discontinued - prias_long$dom_diagnosis)/YEAR_DIVIDER

#Set dom as NA for patient data based on implausible dates of measurements
#2019 june is 13.5 years since PRIAS started in 2006
MAX_FOLLOW_UP_YEARS = 13.5

temp_filter = is.na(prias_long$year_psa) | prias_long$year_psa < 0 | prias_long$year_psa > MAX_FOLLOW_UP_YEARS
prias_long$psa[temp_filter] = NA
prias_long$dom_psa[temp_filter] = NA
prias_long$year_psa[temp_filter] = NA

#same we do for date of DRE and Gleason
temp_filter = is.na(prias_long$year_dre_gleason) | prias_long$year_dre_gleason < 0 | prias_long$year_dre_gleason > MAX_FOLLOW_UP_YEARS
prias_long$dre[temp_filter] = NA
prias_long$dre_recode[temp_filter] = NA
prias_long$dom_dre_gleason[temp_filter] = NA
prias_long$year_dre_gleason[temp_filter] = NA
prias_long$gleason_1[temp_filter] = NA
prias_long$gleason_2[temp_filter] = NA
prias_long$gleason_sum[temp_filter] = NA
prias_long$ncores[temp_filter] = NA
prias_long$ncores_pc[temp_filter] = NA

#Now some PSA are too big to be true, and some are 99, none are negative
summary(prias_long$psa)
sort(prias_long$psa, decreasing = T)[1:50]

temp_filter = is.na(prias_long$psa) | prias_long$psa > 95
prias_long$psa[temp_filter] = NA
prias_long$dom_psa[temp_filter] = NA
prias_long$year_psa[temp_filter] = NA

#Now some PSA are 0.0 and PSA jumps to suddenly zero, such as this guy
View(prias_long[prias_long$P_ID==548, c("psa","year_psa")])

#There are 41 such people
unique(prias_long$P_ID[prias_long$psa %in% 0])

#I checked every row of these 41 people
#PSA is likely really 0 only for patients 1631, 2521, 2565, 2654
#Because first two really had low PSA overall, and last two has 0 PSA after a few years
#I am converting all zero psa to 0.1 as well
View(prias_long[prias_long$P_ID %in% unique(prias_long$P_ID[prias_long$psa %in% 0]), c("P_ID","psa","year_psa")])

temp_filter = prias_long$psa %in% 0 & !(prias_long$P_ID %in% c(1631, 2521, 2565, 2654))
prias_long$psa[temp_filter] = NA
prias_long$dom_psa[temp_filter] = NA
prias_long$year_psa[temp_filter] = NA

#for the 4 aforementioned patients I set PSA 0 as 0.1
prias_long$psa[prias_long$psa %in% 0 & 
                 prias_long$P_ID %in% c(1631, 2521, 2565, 2654)] = 0.1

#Now lets focus on DRE
table(prias_long$dre, useNA = "always")
table(prias_long$dre_recode, useNA = "always")

#So do not use DRE recode as it has too many NA.
#Now remove whitespace from dre levels, and all empty should be set to NA
levels(prias_long$dre) = trimws(levels(prias_long$dre))
prias_long$dre[prias_long$dre %in% c("")] = NA
prias_long$dre = droplevels(prias_long$dre)
table(prias_long$dre, useNA = "always")

#Now lets focus on dummy indicator, dont think of active indicator
temp_filter = prias_long$dummy==1
prias_long$dom_psa[temp_filter] = NA
prias_long$year_psa[temp_filter] = NA
prias_long$psa[temp_filter] = NA

prias_long$dom_dre_gleason[temp_filter] = NA
prias_long$year_dre_gleason[temp_filter] = NA
prias_long$dre[temp_filter] = NA
prias_long$dre_recode[temp_filter] = NA
prias_long$gleason_1[temp_filter] = NA
prias_long$gleason_2[temp_filter] = NA
prias_long$gleason_sum[temp_filter] = NA
prias_long$ncores[temp_filter] = NA
prias_long$ncores_pc[temp_filter] = NA

#Now lets focus on surgery years
summary(prias_long$year_surgery)

table(prias_long$gleason_rad_prost, useNA = "always")

temp_filter = is.na(prias_long$year_surgery) | 
  prias_long$year_surgery<=0 | 
  prias_long$gleason_rad_prost %in% c(NA)

prias_long$Gleason1_Rad_Prost[temp_filter] = NA
prias_long$Gleason2_Rad_Prost[temp_filter] = NA
prias_long$gleason_rad_prost[temp_filter] = NA
prias_long$Surgerydate[temp_filter] = NA
prias_long$year_surgery[temp_filter] = NA

#Now the same for discontinued status and its date
summary(prias_long$year_discontinued)
#None is negative although some are zero, dont worry about them
which(prias_long$year_discontinued %in% 0)

#this is not an issue for dead date
summary(prias_long$year_dead)

#Now lets focus on dates, first for PSA, lets order dates
psa_date_right_order = order(prias_long$P_ID, prias_long$year_psa, na.last = T, decreasing = F)
prias_long$psa = prias_long$psa[psa_date_right_order]
prias_long$year_psa = prias_long$year_psa[psa_date_right_order]
prias_long$dom_psa = prias_long$dom_psa[psa_date_right_order]

#Now ofcourse some people have the psa repeated on same date or almost same date
#I only care if it is measured on the same day twice, so 107 such people are there
psa_date_diff_issue = by(prias_long$year_psa, INDICES = prias_long$P_ID, FUN = function(x){
  any(duplicated(x[!is.na(x)]))
})
table(psa_date_diff_issue)

temp_filter = unlist(by(prias_long$year_psa, INDICES = prias_long$P_ID, FUN = duplicated))
prias_long$psa[temp_filter] = NA
prias_long$year_psa[temp_filter] = NA
prias_long$dom_psa[temp_filter] = NA

#Now lets focus on dates of DRE and gleason, lets order them first
dregleason_date_right_order = order(prias_long$P_ID, prias_long$year_dre_gleason, 
                                    prias_long$gleason_sum, prias_long$dre,
                                    na.last = T,decreasing = F)
prias_long$dom_dre_gleason = prias_long$dom_dre_gleason[dregleason_date_right_order]
prias_long$year_dre_gleason = prias_long$year_dre_gleason[dregleason_date_right_order]
prias_long$dre = prias_long$dre[dregleason_date_right_order]
prias_long$dre_recode = prias_long$dre_recode[dregleason_date_right_order]
prias_long$gleason_1 = prias_long$gleason_1[dregleason_date_right_order]
prias_long$gleason_2 = prias_long$gleason_2[dregleason_date_right_order]
prias_long$gleason_sum = prias_long$gleason_sum[dregleason_date_right_order]
prias_long$ncores = prias_long$ncores[dregleason_date_right_order]
prias_long$ncores_pc = prias_long$ncores_pc[dregleason_date_right_order]

#Now ofcourse some people have the psa repeated on same date or almost same date
#I only care if it is measured on the same day twice, so 355 such people are there
dregleason_date_diff_issue = by(prias_long$year_dre_gleason, 
                                INDICES = prias_long$P_ID, 
                                FUN = function(x){
  sum(duplicated(x[!is.na(x)]))
})
table(dregleason_date_diff_issue)

View(prias_long[prias_long$P_ID %in% as.numeric(names(which(dregleason_date_diff_issue>4))),])

#Remove everything but first duplicated item
temp_filter = unlist(by(prias_long$year_dre_gleason, INDICES = prias_long$P_ID, FUN = duplicated))
prias_long$dom_dre_gleason[temp_filter] = NA
prias_long$year_dre_gleason[temp_filter] = NA
prias_long$dre[temp_filter] = NA
prias_long$dre_recode[temp_filter] = NA
prias_long$gleason_1[temp_filter] = NA
prias_long$gleason_2[temp_filter] = NA
prias_long$gleason_sum[temp_filter] = NA
prias_long$ncores[temp_filter] = NA
prias_long$ncores_pc[temp_filter] = NA

#Now lets focus on gleason reclassification and risk reclassification
#What do we do with Gleason equal to 0?
table(prias_long$gleason_sum)

#In all cases when gleason is zero, number of cores are non zero and non NA
table(prias_long$ncores[prias_long$gleason_sum %in% 0], useNA = 'always')

#But in some cases ncores_pc are more than 0
table(prias_long$ncores_pc[prias_long$gleason_sum %in% 0], useNA = 'always')

#So first lets give Gleason 6 to all those Gleason 0, who have ncores>0 (default) 
#and ncores_pc == 0. that is no cancer found
temp_filter = prias_long$gleason_sum %in% 0 & prias_long$ncores_pc %in% 0
prias_long$gleason_1[temp_filter] = 3
prias_long$gleason_2[temp_filter] = 3
prias_long$gleason_sum[temp_filter] = 6

#now some gleason 0 but non zero ncores_pc. lets check these patients thoroughly
temp_pid_filter=prias_long$P_ID[prias_long$gleason_sum %in% 0 & !(prias_long$ncores %in% c(NA,0))]
View(prias_long[prias_long$P_ID %in% temp_pid_filter,
                c("P_ID","ncores","ncores_pc", "gleason_sum", "year_dre_gleason")])

#After thorough checking only a biopsy for 2520 was an issue. 
#there was not much gap in the biopsies. for the rest set gleason 0 to 6
temp_filter = prias_long$gleason_sum %in% 0 & prias_long$ncores_pc > 0
prias_long$gleason_1[temp_filter & prias_long$P_ID!=2520] = 3
prias_long$gleason_2[temp_filter & prias_long$P_ID!=2520] = 3
prias_long$gleason_sum[temp_filter & prias_long$P_ID!=2520] = 6

prias_long$gleason_1[temp_filter & prias_long$P_ID==2520] = NA
prias_long$gleason_2[temp_filter & prias_long$P_ID==2520] = NA
prias_long$gleason_sum[temp_filter & prias_long$P_ID==2520] = NA

#Now we gotta if two biopsies are too close to each other, 
#because thats not correct, but also doesn't matter coz
#at the end the time interval matters...and that will only be affected
#by a few decimals
biopsy_diff_times = by(prias_long, prias_long$P_ID, FUN = function(x){
  min(diff(x$year_dre_gleason[!is.na(x$gleason_sum)]))
})

biopsy_diff_times_pid = as.numeric(names(sort(biopsy_diff_times)[1:50]))
View(prias_long[prias_long$P_ID %in% biopsy_diff_times_pid & !is.na(prias_long$gleason_sum), 
                c("P_ID", "gleason_sum", "psa","dre","year_dre_gleason")])

#Create new high_gleason variable (takes care of NA)
prias_long$high_gleason = prias_long$gleason_sum > 6

#Now certain patients get a high gleason with a few cores...but they remain in study
#later they undergo biopsy again and they get a low gleason 
#we know Gleason is subject to error. in this case I will assume all Gleason before the
#last high gleason are low gleason, even if the earlier gleasons were low
#there are 18 such patients
sum(by(prias_long$high_gleason, INDICES = prias_long$P_ID, is.unsorted, na.rm=T))

#lets pick up these rows and make them all Gleason 6
temp_filter = unlist(by(prias_long, INDICES = prias_long$P_ID, FUN = function(x){
  #select all except the last one
  if(is.unsorted(x$high_gleason, na.rm = T)){
    return(!is.na(x$high_gleason) & x$year_dre_gleason < max(x$year_dre_gleason[!is.na(x$high_gleason)]))
  }else{
    return(rep(F, 48))
  }
}))

prias_long$gleason_1[temp_filter] = 3
prias_long$gleason_2[temp_filter] = 3
prias_long$gleason_sum[temp_filter] = 6
prias_long$high_gleason[temp_filter] = F

#Now find the Gleason upgrade time for each patient
reclassification_times = by(prias_long, prias_long$P_ID, FUN = function(x){
  high_gleason = x$high_gleason[!is.na(x$high_gleason)]
  year_dre_gleason = x$year_dre_gleason[!is.na(x$high_gleason)]
  
  t_end = min(year_dre_gleason[which(high_gleason==T)])
  # A few patients have a Gleason 6 after they get Gleason 7...clearly measurement error issue
  t_start = max(year_dre_gleason[which(high_gleason==F)])
  
  tright_cens = ifelse(t_end==Inf, t_start, 0.5 * (t_start + t_end))
  year_max_follow_up = ifelse(t_end!=Inf, t_end, 
                              max(x$year_psa[!is.na(x$psa)], 
                                          x$year_dre_gleason[!(is.na(x$dre) & is.na(x$gleason_sum))]))
           
  c(t_start, t_end, tright_cens, year_max_follow_up)
})

reclassification_times = matrix(unlist(reclassification_times), ncol = 4, byrow = T)

prias_long$latest_survival_time = rep(reclassification_times[,1], each=48)
prias_long$earliest_failure_time = rep(reclassification_times[,2], each=48)
prias_long$right_cens_time = rep(reclassification_times[,3], each=48)
prias_long$year_max_followup = rep(reclassification_times[,4], each=48)
prias_long$reclassification = !is.infinite(prias_long$earliest_failure_time)

#some people have follow-up after dying, seems wrong, so fixing it
table(prias_long$year_max_followup > prias_long$year_dead)
prias_long$year_dead[!is.na(prias_long$year_dead) & 
                       prias_long$year_max_followup > prias_long$year_dead] = NA
prias_long$Dead = !is.na(prias_long$year_dead)

#some people have surgery date before maximum follow-up
View(prias_long[!is.na(prias_long$year_surgery) & 
                       prias_long$year_max_followup > prias_long$year_surgery,])

#I checked them in detail, very strange. better leave it the way it is
colnames(prias_long)

SPSS_ORIGIN_DATE = "1582-10-14"
#Now gonna make a single visitTimeYears Variable
prias_long_final = do.call('rbind', by(prias_long, INDICES = prias_long$P_ID, function(x){
  #due to multiple rounds of processing, some NA's remain
  temp_filter = !is.na(x$year_psa) & !is.na(x$psa)
  psa = x$psa[temp_filter]
  dom_psa = x$dom_psa[temp_filter]
  year_psa = x$year_psa[temp_filter]
  
  temp_filter = !is.na(x$year_dre_gleason) & !(is.na(x$dre) & is.na(x$gleason_sum))
  dre = x$dre[temp_filter]
  gleason_sum = x$gleason_sum[temp_filter]
  dom_dre_gleason = x$dom_dre_gleason[temp_filter]
  year_dre_gleason = x$year_dre_gleason[temp_filter]
  
  #Further delete all data after year of maximum follow-up
  year_max_followup = x$year_max_followup[1]
  
  temp_filter = year_psa <= year_max_followup
  dom_psa = dom_psa[temp_filter]
  psa = psa[temp_filter]
  year_psa = year_psa[temp_filter]
  
  temp_filter = year_dre_gleason <= year_max_followup
  dre = dre[temp_filter]
  gleason_sum = gleason_sum[temp_filter]
  dom_dre_gleason = dom_dre_gleason[temp_filter]
  year_dre_gleason = year_dre_gleason[temp_filter]
    
  #creating single year_visit
  new_year_visit = sort(unique(c(year_psa, year_dre_gleason)), decreasing = F)
  new_dom_visit = sort(unique(c(dom_psa, dom_dre_gleason)), decreasing = F)
  
  new_psa = rep(NA, length(new_year_visit))
  new_psa[new_year_visit %in% year_psa] = psa
  
  new_gleason = rep(NA, length(new_year_visit))
  new_dre = factor(rep(NA, length(new_year_visit)), levels=levels(prias_long$dre))
  
  dregleason_indices = new_year_visit %in% year_dre_gleason
  new_gleason[dregleason_indices] = gleason_sum
  new_dre[dregleason_indices] = dre
  
  start_date = as.character.Date(as.POSIXct(x$dom_diagnosis[1],
                                            origin = SPSS_ORIGIN_DATE, "%Y-%m-%d"))
  
  data.frame(P_ID=x$P_ID[1], 
             visit_number=1:length(new_year_visit),
             age=x$age[1],
             dom_diagnosis=x$dom_diagnosis[1],
             start_date=start_date,
             reclassification=x$reclassification[1],
             latest_survival_time = x$latest_survival_time[1],
             earliest_failure_time = x$earliest_failure_time[1],
             right_cens_time = x$right_cens_time[1],
             year_max_followup = x$year_max_followup[1],
             
             Treatment_Type = x$Treatment_Type[1],
             Treatment_Reason = x$Treatment_Reason[1],
             
             year_discontinued = x$year_discontinued[1],
             Dead = x$Dead[1], year_dead = x$year_dead[1], 
             gleason_rad_prost = x$gleason_rad_prost[1], 
             year_surgery = x$year_surgery[1],
             
             year_visit = new_year_visit,
             dom_visit = new_dom_visit,
             psa = new_psa, log2psaplus1 = log(new_psa + 1, base = 2),
             dre = new_dre, palpable_dre = (new_dre!='T1c'),
             gleason_sum = new_gleason, high_gleason = new_gleason > 6)
}))

longitudinal_data = droplevels(prias_long_final)

save(longitudinal_data, file="data/clean_longitudinal_data.Rdata")
