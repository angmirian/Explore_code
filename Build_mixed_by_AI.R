# mixed_by model building for AIC comparison
# Last modified by Angela I Feb 2024
# Mixed by runs an MLM at each time point
# ddf is the output structure; reml and ml are two different convergence methods (use reml - contains AIC and BIC and log likelihood)
# ddf$coeff_df_reml - has all fixed and random effects; value is the estimate and estimated error is std.error

rm(list=ls())

library(stringr)
library(pracma)
library(wesanderson)
library(tidyverse)
library(data.table) 

align_to = "response" #options are clock or response
# AI - these flags allow you to collapse across hemispheres and also to be able to split into networks; default is to do by ROI number (atlas_value)
do_subregion = FALSE #this will collapse across subregions (separate hemispheres)
do_structure = FALSE #collapses across coarser groups and hemispheres
do_network = TRUE #this will collapse across networks and hemispheres
do_hippocampusAP = FALSE #run hippocampus anterior vs. posterior analysis
#models = c(5, 5.1, 5.2, 5.5, 5.6, 5.7, 7, 7.5, 7.6) #models to run, clock
#models = c(5, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 7, 7.5, 7.6) #models to run, feedback
models = c(7,7.5,7.6)
region_to_run = "subcortical" #options are cortical (vmPFC) and subcortical
do_vPFC_HC = FALSE
split_hipp_amyg = TRUE #if want to separate out hippocampus and amygdala

subcortical_cache_dir = '/Users/angela/Documents/RESEARCH/ANALYSES/Medusa_analysis/subcortical'  # where MEDUSA output files are
cortical_cache_dir = '/Users/angela/Documents/RESEARCH/ANALYSES/Medusa_analysis/cortical'
ncores = 3 #I have 8 total on my laptop; put 1 less than total cores on machine or else it will slow down
source("/Users/angela/Documents/RESEARCH/SCRIPTS/fmri.pipeline/R/mixed_by.R")

setwd("~/Documents/RESEARCH/ANALYSES/Explore_behavioral")

# behavioral analyses
#source("~/Documents/RESEARCH/ANALYSES/Explore_behavioral/scripts/get_trial_data.R") #from https://github.com/UNCDEPENdLab/clock_analysis/blob/master/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R
#trial_df <- get_trial_data("~/Documents/RESEARCH/ANALYSES/temporal_instrumental_agent/clock_task/vba_fmri", dataset = "explore") 
#saveRDS(trial_df, file = "trial_df_02142023.rds")
#trial_df <- readRDS("trial_df_02142023.rds")
trial_df <- readRDS("trial_df_03282023.rds") #fixed coding of two subjs
trial_df$id <- as.integer(sub("_1", "", trial_df$id))

#Load subject's data
#sub_df <- readRDS("explore_n146.rds") %>% mutate(id = registration_redcapid)
#pid5_df <- read.csv("~/Documents/RESEARCH/ANALYSES/Explore_behavioral/PID5/pid5_facets.csv")
#sub_df <- left_join(sub_df,pid5_df,by="id")
#demos <- sub_df %>%
#  dplyr::select(-id) %>% # ascending id
#  dplyr::rename(id=registration_redcapid,groupLeth=group_leth,sex=gender) %>% # real study id
#  group_by(id) %>%
#  dplyr::filter(row_number() == 1) %>%
#  ungroup() %>% 
#  dplyr::select(id, age, sex, education_yrs, wtar, exit, mmse, Group, Group_a, groupLeth, age_at_first_attempt, time_since_last_attempt, max_ideation, max_intent, names(pid5_df), starts_with("uppsp"), starts_with("neo"), starts_with("bis")) %>% 
#  droplevels() %>%
#  dplyr::mutate(id = as.integer(id), registration_edu = as.numeric(education_yrs))

sub_df1 <- readRDS("explore_n146.rds") %>% mutate(id = registration_redcapid)
sub_df2 <- read.csv("~/Documents/RESEARCH/ANALYSES/REDCAP_DATA/EXPLORE/explore_clock_indv_diff_new_8_3_23.csv")
sub_df <- left_join(sub_df1,sub_df2,by="registration_redcapid")
pid5_df <- read.csv("~/Documents/RESEARCH/ANALYSES/Explore_behavioral/PID5/pid5_facets.csv")
sub_df <- left_join(sub_df,pid5_df,by="id")
demos <- sub_df %>%
  dplyr::select(-id) %>% # ascending id
  dplyr::rename(id=registration_redcapid,groupLeth=group_leth,sex=gender) %>% # real study id
  group_by(id) %>%
  dplyr::filter(row_number() == 1) %>%
  ungroup() %>% 
  #dplyr::select(id, age, sex, education_yrs, wtar, exit, mmse, Group, Group_a, groupLeth, age_at_first_attempt, time_since_last_attempt, max_ideation, max_intent, names(pid5_df), starts_with("uppsp"), starts_with("neo"), starts_with("bis")) %>% 
  droplevels() %>%
  dplyr::mutate(id = as.integer(id), registration_edu = as.numeric(education_yrs))

#Add in mean score
trial_df <- trial_df %>% group_by(id) %>% mutate(total_score = sum(score_csv)) %>% ungroup()

#recode pid5 to start with "pid5"
demos <- demos %>%
  rename(pid5_anhedonia = "Anhedonia", pid5_anxiousness = "Anxiousness", pid5_attention_seeking = "Attention_Seeking",
         pid5_callousness = "Callousness", pid5_deceitfulness = "Deceitfulness", pid5_depressivity = "Depressivity",
         pid5_distractibility = "Distractibility", pid5_eccentricity = "Eccentricity", pid5_emotional_lability = "Emotional_Lability",
         pid5_grandiosity = "Grandiosity", pid5_hostility = "Hostility", pid5_impulsivity = "Impulsivity", 
         pid5_intimacy_avoidance = "Intimacy_Avoidance", pid5_irresponsibility = "Irresponsibility", pid5_manipulativeness = "Manipulativeness",
         pid5_perceptual_dysregulation = "Perceptual_Dysregulation", pid5_perseveration = "Perseveration", pid5_restricted_affectivity = "Restricted_Affectivity",
         pid5_rigid_perfectionism = "Rigid_Perfectionism", pid5_risk_taking = "Risk_Taking", pid5_separation_insecurity = "Separation_Insecurity",
         pid5_submissiveness = "Submissiveness", pid5_suspiciousness = "Suspiciousness", pid5_unusual_beliefs_experiences = "Unusual_Beliefs_Experiences",
         pid5_withdrawal = "Withdrawal")

####compute subscales:
demos<-demos %>% dplyr::mutate(pid5_negativeaffect = (pid5_anxiousness+ pid5_emotional_lability+ pid5_separation_insecurity)/3)
demos<-demos %>% dplyr::mutate(pid5_detachment = (pid5_withdrawal+ pid5_anhedonia+ pid5_intimacy_avoidance)/3)
demos<-demos %>% dplyr::mutate(pid5_antagonism = (pid5_manipulativeness+ pid5_deceitfulness+ pid5_grandiosity)/3)
demos<-demos %>% dplyr::mutate(pid5_disinhibition = (pid5_irresponsibility+ pid5_impulsivity+ pid5_distractibility)/3)
demos<-demos %>% dplyr::mutate(pid5_psychoticism = (pid5_unusual_beliefs_experiences+ pid5_eccentricity+ pid5_perceptual_dysregulation)/3)
demos<-demos %>% dplyr::mutate(uppsp_total = uppsp_negative_urgency + uppsp_positive_urgency + uppsp_lack_of_premeditation + uppsp_lack_of_perseveration) # AI added 7/25 to get total UPPS score


#Remove excludes 
excludes = c(207224, 210374, 210701, 213708, 216806, 216845, 219392, 220024, 221698, 440149, 440311, 440336) #imaging QC excludes (excluded in original anlaysis)
#excludes = c(207224, 210374, 210701, 213708, 216806, 216845, 219392, 220024, 221698, 440149, 440311, 440336, 211253, 431100, 440243, 440254, 440263, 440392, 221292, 221637, 220913, 431224) #imaging QC excludes, imaging questional QCs, behavioral flatliners, 2 subjs with really long runs
demos <- demos %>% filter(!id %in% excludes)

trial_df <- trial_df %>%
  left_join(demos, by="id") %>%
  mutate(reward_lag_rec = if_else(reward_lag=="omission", -0.5, 0.5),
         reward_lag2 = !omission_lag2,
         reward_lag3 = !omission_lag3,
         reward_lag4 = !omission_lag4)
#remove exclude subjects
trial_df <- trial_df %>% filter(!id %in% excludes)

#Code the condition trial (more similar to bsocial run_trial)
trial_df <- trial_df %>% 
  group_by(id,run_number) %>% 
  arrange(id, run_number, trial) %>% 
  mutate(condition_trial = run_trial-floor(run_trial/40.5)*40,                                                                                          
         condition_trial_neg_inv = -1000 / condition_trial,
         condition_trial_neg_inv_sc = as.vector(scale(condition_trial_neg_inv)),
  ) %>% ungroup()
trial_df$id <- gsub(x = trial_df$id, pattern='881224',replacement='431224')
trial_df$id <- gsub(x = trial_df$id, pattern='881230',replacement='431230')

###########################
#####  AI (1/5/2023)  #####
###########################
message("Loading medusa data from cache: ", subcortical_cache_dir)
if(align_to == "response"){
  print("aligning to response")
  subcortical_df = fread(file.path(subcortical_cache_dir,  'rt_aligned_striatum_hipp_thalamus.csv.gz')) #response aligned
  cortical_df = fread(file.path(cortical_cache_dir,  'rt_aligned_200_vmpfc.csv.gz')) #response aligned
} else if(align_to == "clock"){
  print("aligning to clock")
  subcortical_df = fread(file.path(subcortical_cache_dir,  'clock_aligned_striatum_hipp_thalamus.csv.gz')) #clock aligned
  cortical_df = fread(file.path(cortical_cache_dir,  'clock_aligned_200_vmpfc.csv.gz')) #response aligned
}


#Add vmpfc to the subcortical df
vmpfc_roi_list <- c(55,56,65,66,67,84,86,88,89,159,160,161,170,171,191,192,194)
#vmpfc_roi_list <- append(vmpfc_roi_list,seq(201,230,1))
cortical_df <- cortical_df %>% filter(atlas_value %in% vmpfc_roi_list)

#Merge cortical with subcortical medusa data - takes up too much memory so doesn't work
#subcortical_df <- subcortical_df %>% left_join(cortical_df, by="id") #add vmpfc to subcortical_df 
#df = merge(x = subcortical_df, y = cortical_df, by = "id",
#           all.y = TRUE)

subcortical_df$id <- gsub(x = subcortical_df$id, pattern='881224',replacement='431224')
subcortical_df$id <- gsub(x = subcortical_df$id, pattern='881230',replacement='431230')
cortical_df$id <- gsub(x = cortical_df$id, pattern='881224',replacement='431224')
cortical_df$id <- gsub(x = cortical_df$id, pattern='881230',replacement='431230')

subcortical_df <- subcortical_df %>% rename(sub_decon = decon_mean)
cortical_df <- cortical_df %>% rename(sub_decon = decon_mean)

trial_df <- trial_df %>% 
    group_by(id, run) %>% arrange_by(run_trial) %>%
      mutate(v_max_lead_sc = lead(v_max_wi),
             iti_sc = scale(iti_ideal),
             rt_csv_lag_sc = lag(rt_csv_sc),
             reward_rec = if_else(reward=="omission", -0.5, 0.5),
             reward_lag_rec = if_else(reward_lag=="omission", -0.5, 0.5),
             abs_pe_max_sc = scale(abs(pe_max)),
             pe_max_lag = lag(pe_max),
             pe_max_lag_sc = scale(pe_max_lag),
             abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
             rt_csv_lag_sc = lag(rt_csv_sc),
             iti_prev_sc = lag(iti_sc)) %>% ungroup

#fix how run is coded in subcortical_df
subcortical_df <- subcortical_df %>% mutate(run = case_when(
  str_detect(run, regex("run1", ignore_case=TRUE)) ~ 1,
  str_detect(run, regex("run2", ignore_case=TRUE)) ~ 2))
cortical_df <- cortical_df %>% mutate(run = case_when(
  str_detect(run, regex("run1", ignore_case=TRUE)) ~ 1,
  str_detect(run, regex("run2", ignore_case=TRUE)) ~ 2))
#Add run_trial to subcortical_df
subcortical_df <- subcortical_df %>% mutate(trial = as.vector(trial)) %>% mutate(run_trial = case_when(
trial < 121 ~ trial,
trial > 120 ~ trial-120L)) # the L tells it that I want it to be an integer vector not a double vector

subcortical_df <- subcortical_df %>% filter(evt_time > -3 & evt_time < 3) #cutoff after 5 seconds since this will bleed into other trials
cortical_df <- cortical_df %>% mutate(trial = as.vector(trial)) %>% mutate(run_trial = case_when(
  trial < 121 ~ trial,
  trial > 120 ~ trial-120L)) # the L tells it that I want it to be an integer vector not a double vector
cortical_df <- cortical_df %>% filter(evt_time > -3 & evt_time < 3) #cutoff after 5 seconds since this will bleed into other trials

#Get labels file; change to the adjusted one (region 191 changed from Default to Limbic for symmetry)??
#labels <- read_csv("~/Documents/RESEARCH/ANALYSES/fMRI/Dec2022/extracted_values/12Dec2022/region_labels_244_adj.csv") %>% mutate(roi_num = as.numeric(roi_num)) %>% inner_join(read_csv("~/Documents/RESEARCH/ANALYSES/fMRI/Dec2022/extracted_values/12Dec2022/region_lookup_244.csv"), by = "roi_num")
labels <- read_csv("~/Documents/RESEARCH/ANALYSES/fMRI/Dec2022/extracted_values/12Dec2022/region_labels_244.csv") %>% mutate(roi_num = as.numeric(roi_num)) %>% inner_join(read_csv("~/Documents/RESEARCH/ANALYSES/fMRI/Dec2022/extracted_values/12Dec2022/region_lookup_244.csv"), by = "roi_num")
labels <- labels %>% mutate(atlas_value = as.numeric(roi_num))
# label missing networks as amyg_hipp_thal (amygdala, hippocampus, thalamus)
if (split_hipp_amyg){
  labels$network[str_detect(labels$subregion, regex("Hippocampus", ignore_case=TRUE))] <- "hippocampus"
  labels$network[str_detect(labels$subregion, regex("BLA", ignore_case=TRUE))] <- "amygdala"
  labels$network[str_detect(labels$subregion, regex("CMN", ignore_case=TRUE))] <- "amygdala"
} else {labels$network[is.na(labels$network)] <- "amyg_hipp_thal"}
#Merge labels with subcortical medusa data
subcortical_df <- subcortical_df %>% inner_join(labels, by="atlas_value") #add labels to subcortical_df 
cortical_df <- cortical_df %>% inner_join(labels, by="atlas_value") #add labels to subcortical_df 

#Code the groupings
subcortical_df <- subcortical_df %>% mutate(group = as.factor(case_when(
  str_detect(subregion, regex("Anterior Putamen", ignore_case=TRUE)) ~ "Striatum",
  str_detect(subregion, regex("Caudate Tail and Lateral Putamen", ignore_case=TRUE)) ~ "Striatum",
  str_detect(subregion, regex("Caudate Head", ignore_case=TRUE)) ~ "Striatum",
  str_detect(subregion, regex("Ventral Striatum", ignore_case=TRUE)) ~ "Striatum",
  str_detect(subregion, regex("Posterior Putamen", ignore_case=TRUE)) ~ "Striatum",
  str_detect(subregion, regex("BLA", ignore_case=TRUE)) ~ "Amygdala",
  str_detect(subregion, regex("CMN", ignore_case=TRUE)) ~ "Amygdala",
  str_detect(subregion, regex("Hippocampus", ignore_case=TRUE)) ~ "Hippocampus",
  str_detect(subregion, regex("Thalamus", ignore_case=TRUE)) ~ "Thalamus")))
cortical_df <- cortical_df %>% mutate(group = as.factor(case_when(
  str_detect(subregion, regex("OFC1", ignore_case=TRUE)) ~ "mPFC",
  str_detect(subregion, regex("OFC2", ignore_case=TRUE)) ~ "mPFC",
  str_detect(subregion, regex("OFC3", ignore_case=TRUE)) ~ "mPFC",
  str_detect(subregion, regex("PFC2", ignore_case=TRUE)) ~ "mPFC",
  str_detect(subregion, regex("PFC4", ignore_case=TRUE)) ~ "mPFC",
  str_detect(subregion, regex("PFC6", ignore_case=TRUE)) ~ "mPFC",
  str_detect(subregion, regex("PFC7", ignore_case=TRUE)) ~ "mPFC",
  str_detect(subregion, regex("PFCdPFCm1", ignore_case=TRUE)) ~ "mPFC",
  str_detect(subregion, regex("PFCdPFCm2", ignore_case=TRUE)) ~ "mPFC",
  str_detect(subregion, regex("PFCdPFCm4", ignore_case=TRUE)) ~ "mPFC",
  str_detect(subregion, regex("PFCl1", ignore_case=TRUE)) ~ "mPFC",
  str_detect(subregion, regex("PFCl2", ignore_case=TRUE)) ~ "mPFC")))

#Code hippocampus subregions - anterior vs. posterior
subcortical_df <- subcortical_df %>% mutate(hippocampus_subregion = as.factor(case_when(
  str_detect(subregion, regex("Anterior Putamen", ignore_case=TRUE)) ~ "Striatum",
  str_detect(subregion, regex("Caudate Tail and Lateral Putamen", ignore_case=TRUE)) ~ "Striatum",
  str_detect(subregion, regex("Caudate Head", ignore_case=TRUE)) ~ "Striatum",
  str_detect(subregion, regex("Ventral Striatum", ignore_case=TRUE)) ~ "Striatum",
  str_detect(subregion, regex("Posterior Putamen", ignore_case=TRUE)) ~ "Striatum",
  str_detect(subregion, regex("BLA", ignore_case=TRUE)) ~ "Amygdala",
  str_detect(subregion, regex("CMN", ignore_case=TRUE)) ~ "Amygdala",
  str_detect(subregion, regex("Thalamus", ignore_case=TRUE)) ~ "Thalamus",
  str_detect(subregion, regex("anterior", ignore_case=TRUE)) ~ "AH",
  str_detect(subregion, regex("posterior", ignore_case=TRUE)) ~ "PH")))

#Merge medusa and trial data
if (do_vPFC_HC) {
  Q <- merge(trial_df, cortical_df, by = c("id", "run", "trial","run_trial")) %>% arrange("id","run","run_trial","evt_time") #takes a long time
  Q <- Q %>% rename(vmPFC_decon = decon_mean) %>% select(!decon_median & !decon_sd)
  hc <- subcortical_df %>% filter(atlas_value %in% c(223,224,225,226,227,228,229,230))
  hc <- hc %>% select(!decon_median & !decon_sd)
  hc <- hc %>% mutate(block = case_when(trial <= 40 ~ 1, 
                                        trial > 40 & trial <= 80 ~ 2,
                                        trial > 80 & trial <=120 ~ 3, 
                                        trial > 120 & trial <=160 ~ 4,
                                        trial > 160 & trial <=200 ~ 5,
                                        trial > 200 & trial <=240 ~ 6))
  hc <- hc %>% mutate(HC_region = case_when(atlas_value==223 ~ 'PH',
                                            atlas_value==224 ~ 'PH',
                                            atlas_value==225 ~ 'AH',
                                            atlas_value==226 ~ 'AH',
                                            atlas_value==227 ~ 'PH',
                                            atlas_value==228 ~ 'PH',
                                            atlas_value==229 ~ 'AH',
                                            atlas_value==230 ~ 'AH'))
  hc <- hc %>% select(!atlas_value)
  #Average HC signal
  hc <- hc %>% group_by(id,run,trial,evt_time,HC_region) %>% summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>% ungroup() # 12 -> 2  
  hc <- hc %>% rename(decon_mean=decon1)
  hc <- hc %>% group_by(id,run) %>% mutate(HCwithin = scale(decon_mean),HCbetween=mean(decon_mean,na.rm=TRUE)) %>% ungroup()
  Q <- Q %>% filter(evt_time > -3 & evt_time < 3) #cut from 4 to +/-3 on 6/2/23
  hc <- hc %>% filter(evt_time > -3 & evt_time < 3)
  rm(subcortical_df)
  rm(cortical_df)
  Q <- inner_join(Q,hc,by=c('id','run','trial','run_trial','evt_time'))
  rm(hc)
  
  
} else if(region_to_run == "subcortical") {
  Q <- merge(trial_df, subcortical_df, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time") #takes a long time
  Q <- Q %>% rename(subcortical_decon = sub_decon) %>% filter(atlas_value<231) #accidentally included one cerebellum ROI (231) in mask I sent to Bea
  Q$subcortical_decon[Q$evt_time > Q$iti_ideal] = NA
} else if (region_to_run=="cortical") {
  Q <- merge(trial_df, cortical_df, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time") #takes a long time
  Q <- Q %>% rename(subcortical_decon = sub_decon)
  Q$subcortical_decon[Q$evt_time > Q$iti_ideal] = NA
}

#Censor out data that bleeds into adjacent trials (added 4/21/2023) - I think this is right... but will see when I run the new Medusa models
if(align_to == "response"){
  Q$subcortical_decon[Q$evt_time > Q$iti_ideal] = NA #censors next trial
  Q$subcortical_decon[Q$evt_time < -(Q$rt_csv + Q$iti_prev)] = NA # censors prior trial
} else if(align_to == "clock"){
  Q$subcortical_decon[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA # censors next trial
  Q$subcortical_decon[Q$evt_time < -Q$iti_prev] = NA # censors prior trial (edited) 
}

#To get basic demographics based on lethality group - before scaling age
#Q %>% filter(group_leth_c=="LL_Attempters") %>% select(id,sex,age) %>% unique() %>% summarise(group_n=n(), males_n=sum(sex=="M"),males_percent=100*(males_n/group_n))
#Q %>% filter(group_leth_c=="HL_Attempters") %>% select(id,sex,age) %>% unique() %>% summarise(group_n=n(), males_n=sum(sex=="M"),males_percent=100*(males_n/group_n))


# scale demographics stuff
#demo <- readRDS('/Users/angela/Documents/RESEARCH/ANALYSES/fMRI/Dec2022/extracted_values/12Dec2022/explore_n146.rds') %>% mutate(id = registration_redcapid)
#demos$id <- as.factor(demos$id)
#Q <- inner_join(Q,demos,by=c('id')) # takes a while
Q$female <- relevel(as.factor(Q$sex),ref="F")
Q$age <- scale(Q$age)
Q$wtar <- scale(Q$wtar)
Q$exit <- scale(Q$exit)
Q$education_yrs <- scale(Q$education_yrs)
Q$pid5_psychoticism <- scale(Q$pid5_psychoticism)
Q$uppsp_total <- scale(Q$uppsp_total)
Q$neo_neuroticism <- scale(Q$neo_neuroticism)
Q$drs_total <- scale(Q$drs_total)
Q$cirsg <- scale(Q$cirsg)
Q$ham_hamtotal_17items <- scale(Q$ham_hamtotal_17items)
Q$beckhopelessness <- scale(Q$beckhopelessness)
Q$athf <- scale(Q$athf)
Q <- Q %>% dplyr::rename(group_leth = groupLeth)

#code groupLeth with controls as reference
Q$group_leth_c <- factor(Q$group_leth, level = c("Controls", "Depressed", "Ideators", "LL_Attempters", "HL_Attempters"))
Q$group_leth_LL <- factor(Q$group_leth, level = c("LL_Attempters",  "HL_Attempters", "Controls", "Depressed", "Ideators"))

#Loop through all of the mdoels
for (model_to_run in models) {
  print(model_to_run)
    
  # left hand side of the GLM is outcome, which is specified in the mixed_by variable
  #9/22/2023 - just added extra controls to 5, 5.4, 5.5, 5.9 for complete sensitivity analysis
  #Add variant with all of the control stuff ()
  if(align_to == "response"){ #need to align to v_max_wi_lead if response aligned
    if (model_to_run == 1) {decode_formula <- formula(~ age*v_max_lead_sc + rt_csv_sc + iti_sc + (1|id))} #model 1 - basic model, include RT and ISI to de-noise
    else if (model_to_run == 2) {decode_formula <- formula(~ v_max_lead_sc + rt_csv_sc + rt_csv_lag_sc + iti_sc + reward_rec + reward_lag_rec +  abs_pe_max_sc + (1|id))}
    else if (model_to_run == 3) {decode_formula <- formula(~ v_max_wi + rt_csv_sc + rt_csv_lag_sc + iti_sc + reward_rec + reward_lag_rec +  abs_pe_max_sc + (1|id))} 
    #else if (model_to_run == 2) {decode_formula <- formula(~ rt_csv_sc + rt_csv_lag_sc + iti_sc + reward_rec + reward_lag_rec +  abs_pe_max_sc + (v_max_lead_sc + 1|id))}
    #else if (model_to_run == 3) {decode_formula <- formula(~ rt_csv_sc + rt_csv_lag_sc + iti_sc + reward_rec + reward_lag_rec +  abs_pe_max_sc + (v_max_wi + 1|id))} 
    else if (model_to_run == 4) {decode_formula <- formula(~ Group + v_max_lead_sc + rt_csv_sc + rt_csv_lag_sc + iti_sc + reward_rec + reward_lag_rec +  abs_pe_max_sc + (1|id))} 
    else if (model_to_run == 4.5) {decode_formula <- formula(~ Group_a + v_max_lead_sc + rt_csv_sc + rt_csv_lag_sc + iti_sc + reward_rec + reward_lag_recg +  abs_pe_max_sc + (1|id))}
    else if (model_to_run == 5) {decode_formula <- formula(~ Group*v_max_lead_sc + rt_csv_sc + rt_csv_lag_sc + iti_sc + reward_rec + reward_lag_rec +  abs_pe_max_sc + iti_prev_sc + age + (1|id))} #added from + iti_prev_sc on 
    #else if (model_to_run == 5) {decode_formula <- formula(~ Group*v_max_lead_sc + rt_csv_sc + rt_csv_lag_sc + iti_sc + reward_rec + reward_lag_rec +  abs_pe_max_sc + iti_prev_sc + age + wtar + education_yrs + drs_total + cirsg + ham_hamtotal_17items + beckhopelessness + athf + anxiety_dx + (1|id))} #added from + iti_prev_sc on 
    else if (model_to_run == 5.1) {decode_formula <- formula(~ Group*v_max_lead_sc + rt_csv_sc + rt_csv_lag_sc + iti_sc + reward_rec + reward_lag_rec +  abs_pe_max_sc + iti_prev_sc + age + (1|id))} # add age and iti_prev
    else if (model_to_run == 5.2) {decode_formula <- formula(~ Group*v_max_lead_sc + rt_csv_sc + rt_csv_lag_sc + iti_sc + reward_rec + reward_lag_rec +  abs_pe_max_sc + wtar + (1|id))} 
    else if (model_to_run == 5.3) {decode_formula <- formula(~ Group*v_max_wi + rt_csv_sc + rt_csv_lag_sc + iti_sc + reward_rec + reward_lag_rec +  abs_pe_max_sc + (1|id))} 
    else if (model_to_run == 5.4) {decode_formula <- formula(~ Group*v_max_lead_sc + Group*v_max_wi + rt_csv_sc + rt_csv_lag_sc + iti_sc + reward_rec + reward_lag_rec +  abs_pe_max_sc + iti_prev_sc + age + wtar + education_yrs + drs_total + cirsg + ham_hamtotal_17items + beckhopelessness + athf + anxiety_dx + (1|id))} #both prev and updated value in model
    else if (model_to_run == 5.5) {decode_formula <- formula(~ Group_a*v_max_lead_sc + rt_csv_sc + rt_csv_lag_sc + iti_sc + reward_rec + reward_lag_rec +  abs_pe_max_sc + iti_prev_sc + age + (1|id))} 
    #else if (model_to_run == 5.5) {decode_formula <- formula(~ Group_a*v_max_lead_sc + rt_csv_sc + rt_csv_lag_sc + iti_sc + reward_rec + reward_lag_rec +  abs_pe_max_sc + iti_prev_sc + age + wtar + education_yrs + drs_total + cirsg + ham_hamtotal_17items + beckhopelessness + athf + anxiety_dx + (1|id))} 
    else if (model_to_run == 5.6) {decode_formula <- formula(~ Group_a*v_max_lead_sc + rt_csv_sc + rt_csv_lag_sc + iti_sc + reward_rec + reward_lag_rec +  abs_pe_max_sc + iti_prev_sc + age + (1|id))} #added age and iti_prev_sc
    else if (model_to_run == 5.7) {decode_formula <- formula(~ Group_a*v_max_lead_sc + rt_csv_sc + rt_csv_lag_sc + iti_sc + reward_rec + reward_lag_rec +  abs_pe_max_sc + wtar + (1|id))} 
    else if (model_to_run == 5.8) {decode_formula <- formula(~ Group_a*v_max_wi + rt_csv_sc + rt_csv_lag_sc + iti_sc + reward_rec + reward_lag_rec +  abs_pe_max_sc + (1|id))} # model 5.5 using previous trial's value
    else if (model_to_run == 5.9) {decode_formula <- formula(~ Group_a*v_max_lead_sc + Group_a*v_max_wi + rt_csv_sc + rt_csv_lag_sc + iti_sc + reward_rec + reward_lag_rec +  abs_pe_max_sc + iti_prev_sc + age + wtar + education_yrs + drs_total + cirsg + ham_hamtotal_17items + beckhopelessness + athf + anxiety_dx + (1|id))} # model 5.5 using previous trial's value
    else if (model_to_run == 6) {decode_formula <- formula(~ wtar*age*v_max_lead_sc + rt_csv_sc + iti_sc + (1|id))} #model 1 with wtar interaction
    else if (model_to_run == 6.5) {decode_formula <- formula(~ exit*age*v_max_lead_sc + rt_csv_sc + iti_sc + (1|id))} 
    else if (model_to_run == 7) {decode_formula <- formula(~ group_leth_c*v_max_lead_sc + rt_csv_sc + rt_csv_lag_sc + iti_sc + reward_rec + reward_lag_rec +  abs_pe_max_sc + (1|id))} 
    else if (model_to_run == 7.5) {decode_formula <- formula(~ group_leth*v_max_lead_sc + rt_csv_sc + rt_csv_lag_sc + iti_sc + reward_rec + reward_lag_rec +  abs_pe_max_sc + (1|id))} 
    else if (model_to_run == 7.6) {decode_formula <- formula(~ group_leth_LL*v_max_lead_sc + rt_csv_sc + rt_csv_lag_sc + iti_sc + reward_rec + reward_lag_rec +  abs_pe_max_sc + (1|id))} 
    # Add model with rt_swings behavior if this is possible
    #  else if (model_to_run == 8) {decode_formula <- formula(~ Group*v_max_lead_sc + rt_csv_sc + rt_csv_lag_sc + iti_sc + reward_rec + reward_lag_rec +  abs_pe_max_sc + (1|id))} 
  } else if(align_to == "clock"){
    if (model_to_run == 1) {decode_formula <- formula(~ age*v_max_wi + rt_csv_sc + iti_sc + (1|id))} #model 1 - basic model, include RT and ISI to de-noise
    else if (model_to_run == 1.5) {decode_formula <- formula(~ age*v_max_wi + rt_csv_lag_sc + iti_prev_sc + (1|id))} #- use lagged rtcsv and iti variables
    #else if (model_to_run == 2) {decode_formula <- formula(~ v_max_wi + rt_csv_sc + rt_csv_lag_sc + iti_sc + reward_rec + reward_lag_rec +  abs_pe_max_sc + (1|id))}
    else if (model_to_run == 2) {decode_formula <- formula(~ v_max_wi + rt_csv_lag_sc + iti_prev_sc + reward_lag_rec +  abs_pe_max_lag_sc + (1|id))}
    else if (model_to_run == 3) {stop('REDUNDANT - THIS IS THE SAME AS MODEL 2 for CLOCK-ALIGNED DATA')}
    else if (model_to_run == 4) {decode_formula <- formula(~ Group + v_max_wi + rt_csv_lag_sc + iti_prev_sc + reward_lag_rec +  abs_pe_max_lag_s +  abs_pe_max_sc + (1|id))} 
    else if (model_to_run == 4.5) {decode_formula <- formula(~ Group_a + v_max_wi + rt_csv_lag_sc + iti_prev_sc + reward_lag_rec +  abs_pe_max_lag_sc + (1|id))} 
    else if (model_to_run == 5) {decode_formula <- formula(~ Group*v_max_wi + rt_csv_lag_sc + iti_prev_sc + reward_lag_rec +  abs_pe_max_lag_sc + iti_prev_sc + age + (1|id))} 
    #else if (model_to_run == 5) {decode_formula <- formula(~ Group*v_max_wi + rt_csv_lag_sc + iti_prev_sc + reward_lag_rec +  abs_pe_max_lag_sc + iti_prev_sc + age + wtar + education_yrs + drs_total + cirsg + ham_hamtotal_17items + beckhopelessness + athf + anxiety_dx + (1|id))} 
    else if (model_to_run == 5.1) {decode_formula <- formula(~ Group*v_max_wi + rt_csv_lag_sc + iti_prev_sc + reward_lag_rec +  abs_pe_max_lag_sc + age + (1|id))} # add age to model 5
    else if (model_to_run == 5.2) {decode_formula <- formula(~ Group*v_max_wi + rt_csv_lag_sc + iti_prev_sc + reward_lag_rec +  abs_pe_max_lag_sc + wtar + (1|id))} 
    else if (model_to_run == 5.5) {decode_formula <- formula(~ Group_a*v_max_wi + rt_csv_lag_sc + iti_prev_sc + reward_lag_rec +  abs_pe_max_lag_sc + iti_prev_sc + age + (1|id))} 
    #else if (model_to_run == 5.5) {decode_formula <- formula(~ Group_a*v_max_wi + rt_csv_lag_sc + iti_prev_sc + reward_lag_rec +  abs_pe_max_lag_sc + iti_prev_sc + age + wtar + education_yrs + drs_total + cirsg + ham_hamtotal_17items + beckhopelessness + athf + anxiety_dx + (1|id))} 
    else if (model_to_run == 5.6) {decode_formula <- formula(~ Group_a*v_max_wi + rt_csv_lag_sc + iti_prev_sc + reward_lag_rec +  abs_pe_max_lag_sc + age + (1|id))} #add iti_prev_sc and age
    else if (model_to_run == 5.7) {decode_formula <- formula(~ Group_a*v_max_wi + rt_csv_lag_sc + iti_prev_sc + reward_lag_rec +  abs_pe_max_lag_sc + wtar + (1|id))} #add wtar
    else if (model_to_run == 6) {decode_formula <- formula(~ wtar*age*v_max_wi + rt_csv_lag_sc + iti_prev_sc + (1|id))} #- use lagged rtcsv and iti variables
    else if (model_to_run == 6.5) {decode_formula <- formula(~ exit*age*v_max_wi + rt_csv_lag_sc + iti_prev_sc + (1|id))} #- use lagged rtcsv and iti variables
    else if (model_to_run == 7) {decode_formula <- formula(~ group_leth_c*v_max_wi + rt_csv_lag_sc + iti_prev_sc + reward_lag_rec +  abs_pe_max_lag_sc + (1|id))} 
    else if (model_to_run == 7.5) {decode_formula <- formula(~ group_leth*v_max_wi + rt_csv_lag_sc + iti_prev_sc + reward_lag_rec +  abs_pe_max_lag_sc + (1|id))} 
    else if (model_to_run == 7.6) {decode_formula <- formula(~ group_leth_LL*v_max_wi + rt_csv_lag_sc + iti_prev_sc + reward_lag_rec +  abs_pe_max_lag_sc + (1|id))} 
  }
  
  #Andrew's models
  #decode_formula[[1]] = formula(~ age (1|id/run))
  #decode_formula[[1]] = formula(~ age + female + v_entropy_sc + v_entropy_wi_change + outcome + trial_neg_inv_sc + v_max_wi + rt_csv_sc + iti_sc + rt_vmax_change_sc + (1|id/run))
  #decode_formula[[2]] = formula(~ age + female + v_entropy_sc + v_entropy_wi_change + outcome + trial_neg_inv_sc + v_max_wi + rt_csv_sc + iti_sc + rt_vmax_change_sc + (1 + v_entropy_sc |id/run))
  #decode_formula[[2]] = formula(~ v_entropy_sc*last_outcome + v_entropy_sc*trial_neg_inv_sc + v_max_wi*last_outcome + v_entropy_wi_change + rt_csv_sc + iti_sc + (1|id/run))
  #decode_formula[[2]] = formula(~ v_entropy_sc*last_outcome + v_entropy_sc*trial_neg_inv_sc + v_max_wi*last_outcome + v_entropy_wi_change + rt_csv_sc + iti_sc +  (1+v_max_wi + v_entropy_sc |id/run))
  #if (do_symmetry){ #does L and R of each region together, assuming they don't differ
  
  if (do_structure){
    splits = c('evt_time','group') #will do a separate MLM for each group in splits; evt_time is the TR; collapse across both hemispheres of subregions
    if(region_to_run == "subcortical") {dir <- paste0(subcortical_cache_dir,'/Structure_CombinedHemis')
    dir.create(file.path(subcortical_cache_dir, '/Structure_CombinedHemis'), showWarnings = FALSE)
    } else if(region_to_run == "cortical") {dir <- paste0(cortical_cache_dir,'/Structure_CombinedHemis')
    dir.create(file.path(cortical_cache_dir, '/Structure_CombinedHemis'), showWarnings = FALSE)}
    setwd(dir)
    df0 <- decode_formula   #df0 <- decode_formula[[i]]
    print(df0)
    ddf <- mixed_by(Q, outcomes = "subcortical_decon", rhs_model_formulae = df0 , split_on = splits,
                    padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3, #ncores put 1 less than total cores on your machine or else it will slow down
                    tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE), #specify what output dataframe will have; will be looking mostly at fixed terms; lmer summary, will need all for that
                    #Run post-hoc tests
                    emmeans_spec = list(
                      if(align_to == "response"){ #need to align to v_max_lead_sc if response aligned
                        if (model_to_run == 3 | model_to_run == 5.3 | model_to_run == 5.8) {V = list(outcome="subcortical_decon", model_name="model1",specs=c("v_max_wi"), at=list(v_max_wi=c(-1.5,1.5)))}
                        else {V = list(outcome="subcortical_decon", model_name="model1",specs=c("v_max_lead_sc"), at=list(v_max_lead_sc=c(-1.5,1.5)))}
                      } else if(align_to == "clock"){
                        V = list(outcome="subcortical_decon", model_name="model1",
                                 specs=c("v_max_wi"), at=list(v_max_wi=c(-1.5,1.5)))
                        }
                      )
                    )
  } else if (do_subregion){
    splits = c('evt_time','subregion') #will do a separate MLM for each group in splits; evt_time is the TR; subregions - collapsed across hemispheres
    dir <- paste0(subcortical_cache_dir,'/Subregion_CombinedHemis')
    dir.create(file.path(subcortical_cache_dir, '/Subregion_CombinedHemis'), showWarnings = FALSE)
    setwd(dir)
    df0 <- decode_formula   #df0 <- decode_formula[[i]]
    print(df0)
    ddf <- mixed_by(Q, outcomes = "subcortical_decon", rhs_model_formulae = df0 , split_on = splits,
            padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3, #ncores put 1 less than total cores on your machine or else it will slow down
            tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE), #specify what output dataframe will have; will be looking mostly at fixed terms; lmer summary, will need all for that
            #Run post-hoc tests
            emmeans_spec = list(
                if(align_to == "response"){ #need to align to v_max_lead_sc if response aligned
                  if (model_to_run == 3 | model_to_run == 5.3 | model_to_run == 5.8) {V = list(outcome="subcortical_decon", model_name="model1",specs=c("v_max_lead_sc"), at=list(v_max_lead_sc=c(-1.5,1.5)))}
                  else {V = list(outcome="subcortical_decon", model_name="model1",specs=c("v_max_lead_sc"), at=list(v_max_lead_sc=c(-1.5,1.5)))}
                } else if(align_to == "clock"){
                    V = list(outcome="subcortical_decon", model_name="model1",specs=c("v_max_wi"), at=list(v_max_wi=c(-1.5,1.5)))}
                )
            )
  } else if (do_network){
    splits = c('evt_time','network') #will do a separate MLM for each group in splits; evt_time is the TR; collapse networks; but several are NA (change to amygdala_thalamus_hippocampus)
    if(region_to_run == "subcortical") {dir <- paste0(subcortical_cache_dir,'/Network_CombinedHemis')
      dir.create(file.path(subcortical_cache_dir, '/Network_CombinedHemis'), showWarnings = FALSE)
      #AI added 5/16 to narrow down networks for subcortical models to just DMN, CONT, LIM, and Hipp/Amyg
      if (split_hipp_amyg){Q <- Q %>% filter(network %in% c("amygdala","hippocampus", "Cont","Default","Limbic"))
      } else {Q <- Q %>% filter(network %in% c("amyg_hipp_thal","Cont","Default","Limbic")) %>% filter(!group %in% c("Thalamus"))}
    } else if(region_to_run == "cortical") {dir <- paste0(cortical_cache_dir,'/Network_CombinedHemis')
      dir.create(file.path(cortical_cache_dir, '/Network_CombinedHemis'), showWarnings = FALSE)}
    setwd(dir)
    df0 <- decode_formula   #df0 <- decode_formula[[i]]
    print(df0)
    ddf <- mixed_by(Q, outcomes = "subcortical_decon", rhs_model_formulae = df0 , split_on = splits,
            padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3, #ncores put 1 less than total cores on your machine or else it will slow down
            tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE), #specify what output dataframe will have; will be looking mostly at fixed terms; lmer summary, will need all for that
            #Run post-hoc tests
            emmeans_spec = list(
              if(align_to == "response"){ #need to align to v_max_lead_sc if response aligned
                if (model_to_run == 3 | model_to_run == 5.3 | model_to_run == 5.8) {V = list(outcome="subcortical_decon", model_name="model1",specs=c("v_max_lead_sc"), at=list(v_max_lead_sc=c(-1.5,1.5)))}
                else {V = list(outcome="subcortical_decon", model_name="model1",specs=c("v_max_lead_sc"), at=list(v_max_lead_sc=c(-1.5,1.5)))}
              } else if(align_to == "clock"){
                V = list(outcome="subcortical_decon", model_name="model1",
                         specs=c("v_max_wi"), at=list(v_max_wi=c(-1.5,1.5)))
                }
              )
            )
  } else if (do_hippocampusAP) {
    splits = c('evt_time','hippocampus_subregion') #will do a separate MLM for each group in splits; evt_time is the TR; subregions - collapsed across hemispheres
    dir <- paste0(subcortical_cache_dir,'/Structure_splitHippocampus_CombinedHemis')
    dir.create(file.path(subcortical_cache_dir, '/Structure_splitHippocampus_CombinedHemis'), showWarnings = FALSE)
    setwd(dir)
    df0 <- decode_formula 
    print(df0)
    ddf <- mixed_by(Q, outcomes = "subcortical_decon", rhs_model_formulae = df0 , split_on = splits,
                    padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3, #ncores put 1 less than total cores on your machine or else it will slow down
                    tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE), #specify what output dataframe will have; will be looking mostly at fixed terms; lmer summary, will need all for that
                    #Run post-hoc tests
                    emmeans_spec = list(
                      if(align_to == "response"){ #need to align to v_max_lead_sc if response aligned
                        if (model_to_run == 3 | model_to_run == 5.3 | model_to_run == 5.8) {V = list(outcome="subcortical_decon", model_name="model1",specs=c("v_max_wi"), at=list(v_max_lead_sc=c(-1.5,1.5)))}
                        else {V = list(outcome="subcortical_decon", model_name="model1",specs=c("v_max_lead_sc"), at=list(v_max_lead_sc=c(-1.5,1.5)))}
                      } else if(align_to == "clock"){
                        V = list(outcome="subcortical_decon", model_name="model1",specs=c("v_max_wi"), at=list(v_max_wi=c(-1.5,1.5)))}
                    )
    )
  } else {
    splits = c('evt_time','atlas_value') #will do a separate MLM for each group in splits; evt_time is the TR
    dir <- paste0(subcortical_cache_dir,'/AtlasValue_SeparateHemis')
    dir.create(file.path(subcortical_cache_dir, '/AtlasValue_SeparateHemis'), showWarnings = FALSE)
    setwd(dir)
    df0 <- decode_formula   #df0 <- decode_formula[[i]]
    print(df0)
      ddf <- mixed_by(Q, outcomes = "subcortical_decon", rhs_model_formulae = df0 , split_on = splits,
            padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3, #ncores put 1 less than total cores on your machine or else it will slow down
            tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE), #specify what output dataframe will have; will be looking mostly at fixed terms; lmer summary, will need all for that
            #Run post-hoc tests
            emmeans_spec = list(
              if(align_to == "response"){ #need to align to v_max_lead_sc if response aligned
                if (model_to_run == 3 | model_to_run == 5.3 | model_to_run == 5.8) {V = list(outcome="subcortical_decon", model_name="model1",specs=c("v_max_wi"), at=list(v_max_lead_sc=c(-1.5,1.5)))}
                else {V = list(outcome="subcortical_decon", model_name="model1",specs=c("v_max_lead_sc"), at=list(v_max_lead_sc=c(-1.5,1.5)))}
              } else if(align_to == "clock"){
                V = list(outcome="subcortical_decon", model_name="model1",
                         specs=c("v_max_wi"), at=list(v_max_wi=c(-1.5,1.5)))
                }
              )
            )
  }
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  #need to code in way to change the name of this file based on the model
  if(align_to == "response"){
    if(region_to_run == "subcortical") {save(ddf,file=paste0(curr_date,'-subcortical-response-','model_',model_to_run,'.Rdata')) # removed i (i is for if you are running several models)
    } else if (region_to_run == "cortical") {save(ddf,file=paste0(curr_date,'-vmpfc-response-','model_',model_to_run,'.Rdata'))} # removed i (i is for if you are running several models)
  } else if(align_to == "clock"){
    if(region_to_run == "subcortical") {save(ddf,file=paste0(curr_date,'-subcortical-clock-','model_',model_to_run,'.Rdata')) # removed i (i is for if you are running several models)
    } else if (region_to_run == "cortical") {save(ddf,file=paste0(curr_date,'-vmpfc-clock-','model_',model_to_run,'.Rdata'))} # removed i (i is for if you are running several models)
  }
}

#Now run plot_mixed_by_AI.R
