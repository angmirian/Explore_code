# Modular code to run mixed_by
# 26 April 2024
# Angela Ianni

# Mixed by runs an MLM at each time point
# ddf is the output structure; reml and ml are two different convergence methods (use reml - contains AIC and BIC and log likelihood)
# ddf$coeff_df_reml - has all fixed and random effects; value is the estimate and estimated error is std.error

rm(list=ls())

library(stringr)
library(pracma)
library(wesanderson)
library(tidyverse)
library(data.table) 

########## User Settings ##########
align_to = "clock" #options are clock or response
structure = "network"
region_to_run = "cortical" #options are cortical (vmPFC) and subcortical
subjs_to_run = "all" #options are all, controls_only, or patients_only
trial_data = ("/Users/angela/Documents/Research/Explore/Medusa_analysis/trial_df_03282023.rds")
subject_demographics = ("/Users/angela/Documents/Research/Explore/Medusa_analysis/explore_complete_demographics.rds")
excludes = c(207224, 210374, 210701, 213708, 216806, 216845, 219392, 220024, 221698, 440149, 440311, 440336) #imaging QC excludes (excluded in original anlaysis); excludes = c(207224, 210374, 210701, 213708, 216806, 216845, 219392, 220024, 221698, 440149, 440311, 440336, 211253, 431100, 440243, 440254, 440263, 440392, 221292, 221637, 220913, 431224) #imaging QC excludes, imaging questional QCs, behavioral flatliners, 2 subjs with really long runs
select_networks=TRUE #can only set to TRUE for subcortical; if true only do Cont, Default, Limbic, hipp, and amygdala
if (region_to_run=="cortical") {
  select_networks=FALSE
}

filter_first_10_trials=FALSE #to match what Andrew did - NOPE, NOT WHAT HE DID; IGNORE THIS

ncores = 3 #I have 8 total on my laptop; put 1 less than total cores on machine or else it will slow down
source("/Users/angela/Documents/Research/Scripts/fmri.pipeline/R/mixed_by.R")
setwd("/Users/angela/Documents/Research/Explore/Medusa_analysis")
#Define how output files are named
if (subjs_to_run=="all") {
  if (select_networks==TRUE) {
    output_file_specifier = "_SelectNetworks_RANDOMSLOPES"
  } else {output_file_specifier = "_RANDOMSLOPES"}
} else if (subjs_to_run=="controls_only") {
  if (select_networks==TRUE) {
    output_file_specifier = "_CONTROLS_ONLY_SelectNetworks_RANDOMSLOPES"
  } else {output_file_specifier = "_CONTROLS_ONLY_RANDOMSLOPES"}
} else if (subjs_to_run=="patients_only") {
  if (select_networks==TRUE) {
    output_file_specifier = "_PATIENTS_ONLY_SelectNetworks_RANDOMSLOPES"
  } else {output_file_specifier = "_PATIENTS_ONLY_RANDOMSLOPES"}
}
#Define subset of data to networks to run
if (region_to_run=="cortical" && select_networks==TRUE) {select_networks=FALSE
} else if (region_to_run=="subcortical" && select_networks==TRUE) {
  network_subset=c("Amygdala","Cont","Default","Limbic","Hippocampus_cobra")
}

########## Read in models to run ##########
if (align_to=="clock") {
  if (subjs_to_run=="all") {
    models<-read.table("/Users/angela/Documents/Research/Explore/Medusa_analysis/Decode_formulas/clock_decode_formulas_30May2024_impulsivity_WSLS.txt")
    #models<-read.table("/Users/angela/Documents/Research/Explore/Medusa_analysis/Decode_formulas/clock_decode_formulas_1May2024.txt")
  } else if (subjs_to_run=="controls_only") {
    models<-read.table("/Users/angela/Documents/Research/Explore/Medusa_analysis/Decode_formulas/clock_decode_formulas_CONTROLS_ONLY_3May2024.txt")
  } else if (subjs_to_run=="patients_only") {
    models<-read.table("/Users/angela/Documents/Research/Explore/Medusa_analysis/Decode_formulas/clock_decode_formulas_30May2024_impulsivity_WSLS.txt")
  }

} else if (align_to=="response") {
  if (subjs_to_run=="all") {
    models<-read.table("/Users/angela/Documents/Research/Explore/Medusa_analysis/Decode_formulas/feedback_decode_formulas_30May2024_impulsivity_WSLS.txt")
    #models<-read.table("/Users/angela/Documents/Research/Explore/Medusa_analysis/Decode_formulas/feedback_decode_formulas_1May2024.txt")
  } else if (subjs_to_run=="controls_only") {
    models<-read.table("/Users/angela/Documents/Research/Explore/Medusa_analysis/Decode_formulas/feedback_decode_formulas_CONTROLS_ONLY_3May2024.txt")
  } else if (subjs_to_run=="patients_only") {
    models<-read.table("/Users/angela/Documents/Research/Explore/Medusa_analysis/Decode_formulas/feedback_decode_formulas_30May2024_impulsivity_WSLS.txt")
  }
}
decode_formula = list()
for (i in 1:length(models$V1)) {
  decode_formula[[i]] <- formula(models$V1[[i]])
}
  


########## Load trial and demographic data ##########
trial_df <- readRDS(trial_data)
trial_df$id <- as.integer(sub("_1", "", trial_df$id))
load('/Users/angela/Documents/Research/Explore/Medusa_analysis/combined_rs_betas.Rda')
# Recode omissions as -0.5 and 0.5 and code the condition trial (more similar to bsocial run_trial)
trial_df <- trial_df %>% 
  group_by(id,run_number) %>% 
  arrange(id, run_number, trial) %>% 
  mutate(reward_lag_rec = if_else(reward_lag=="omission", -0.5, 0.5),
         condition_trial = run_trial-floor(run_trial/40.5)*40,                                                                                          
         condition_trial_neg_inv = -1000 / condition_trial,
         condition_trial_neg_inv_sc = as.vector(scale(condition_trial_neg_inv)),
  ) %>% ungroup()
demos <- readRDS(subject_demographics)
#Merge demos and trial_df
trial_df <- trial_df %>% left_join(demos, by="id")
#Merge WSLS behavior and trial_df
wsls_df <- combined_rs_betas %>% select(id,WSLS_rs,LS_rs)
wsls_df$id <- as.integer(wsls_df$id)
trial_df <- trial_df %>% left_join(wsls_df, by="id")
#Remove excluded subjects
trial_df <- trial_df %>% filter(!id %in% excludes)
#Fix mis-named IDs for Explore
trial_df$id <- gsub(x = trial_df$id, pattern='881224',replacement='431224')
trial_df$id <- gsub(x = trial_df$id, pattern='881230',replacement='431230')
#Add extra parameters to trial_df
trial_df <- trial_df %>% 
  group_by(id, run) %>% arrange(id, run, trial) %>% 
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



########## Load Schaeffer Labels ########## 
#Note: adjusted: (region 191 changed from Default to Limbic for symmetry)
labels <- read_csv("/Users/angela/Documents/Research/Explore/fMRI/Dec2022/extracted_values/12Dec2022/region_labels_244.csv") %>% mutate(roi_num = as.numeric(roi_num)) %>% inner_join(read_csv("/Users/angela/Documents/Research/Explore/fMRI/Dec2022/extracted_values/12Dec2022/region_lookup_244.csv"), by = "roi_num")
labels <- labels %>% mutate(atlas_value = as.numeric(roi_num))

# label missing networks as amygdala, hippocampus, thalamus
labels$network[str_detect(labels$subregion, regex("Hippocampus", ignore_case=TRUE))] <- "Hippocampus"
labels$network[str_detect(labels$subregion, regex("BLA", ignore_case=TRUE))] <- "Amygdala"
labels$network[str_detect(labels$subregion, regex("CMN", ignore_case=TRUE))] <- "Amygdala"
labels$network[str_detect(labels$subregion, regex("Thalamus", ignore_case=TRUE))] <- "Thalamus"



########## Load Medusa Data ########## 
#Define data location
subcortical_cache_dir="/Users/angela/Documents/Research/Explore/Medusa_analysis/subcortical"
cortical_cache_dir="/Users/angela/Documents/Research/Explore/Medusa_analysis/cortical"
if (align_to=="clock") {
  print("aligning to clock")
  HCmedusa_data = c("/Users/angela/Documents/Research/Explore/Medusa_analysis/Medusa_data/12slice_cobra_hippocampus_mask/Explore_HC_l_clock.csv.gz",
                  "/Users/angela/Documents/Research/Explore/Medusa_analysis/Medusa_data/12slice_cobra_hippocampus_mask/Explore_HC_r_clock.csv.gz")
  subcortical_data = "/Users/angela/Documents/Research/Explore/Medusa_analysis/subcortical/clock_aligned_striatum_hipp_thalamus.csv.gz"
  cortical_data = "/Users/angela/Documents/Research/Explore/Medusa_analysis/cortical/clock_aligned_200_vmpfc.csv.gz"
} else if (align_to=="response") {
  print("aligning to response")
  HCmedusa_data = c("/Users/angela/Documents/Research/Explore/Medusa_analysis/Medusa_data/12slice_cobra_hippocampus_mask/Explore_HC_l_fb.csv.gz",
                  "/Users/angela/Documents/Research/Explore/Medusa_analysis/Medusa_data/12slice_cobra_hippocampus_mask/Explore_HC_r_fb.csv.gz")
  subcortical_data = "/Users/angela/Documents/Research/Explore/Medusa_analysis/subcortical/rt_aligned_striatum_hipp_thalamus.csv.gz"
  cortical_data = "/Users/angela/Documents/Research/Explore/Medusa_analysis/cortical/rt_aligned_200_vmpfc.csv.gz"
}
#Load the data
if (region_to_run=="cortical") {
  cortical_df = fread(cortical_data)
  cortical_df$id <- gsub(x = cortical_df$id, pattern='881224',replacement='431224')
  cortical_df$id <- gsub(x = cortical_df$id, pattern='881230',replacement='431230')
  vmpfc_roi_list <- c(55,56,65,66,67,84,86,88,89,159,160,161,170,171,191,192,194)
  cortical_df <- cortical_df %>% filter(atlas_value %in% vmpfc_roi_list)
  cortical_df <- cortical_df
  #fix how run is coded
  cortical_df <- cortical_df %>% mutate(run = case_when(
    str_detect(run, regex("run1", ignore_case=TRUE)) ~ 1,
    str_detect(run, regex("run2", ignore_case=TRUE)) ~ 2))
  #Add run_trial to medusa_df
  cortical_df <- cortical_df %>% mutate(trial = as.vector(trial)) %>% mutate(run_trial = case_when(
    trial < 121 ~ trial,
    trial > 120 ~ trial-120L)) # the L tells it that I want it to be an integer vector not a double vector
  #Filter data +/- 3 sec from event
  cortical_df <- cortical_df %>% filter(evt_time > -3 & evt_time < 3)
  #Merge labels with medusa data
  cortical_df <- cortical_df %>% inner_join(labels, by="atlas_value")
  #Code the region groupings
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
  #Merge trial_df and medusa data
  Q <- merge(trial_df, cortical_df, by = c("id", "run", "run_trial")) %>% dplyr::arrange("id","run","run_trial","evt_time") #takes a long time
  Q$decon_mean[Q$evt_time > Q$iti_ideal] = NA
} else if (region_to_run=="subcortical") {
  HCdata1=fread(file.path(HCmedusa_data[1])) # Left; atlas values for HCdata1 and 2 are between 0 and 1
  HCdata2=fread(file.path(HCmedusa_data[2])) # Right
  subcortical_df = fread(subcortical_data)
  #Average across hippocampus slices for each hemisphere
  avg_Lhipp <- HCdata1 %>% group_by(id,run,trial,evt_time) %>% 
    dplyr::summarize(Lhipp_avg_decon_mean=mean(decon_mean,na.rm=TRUE),
                     Lhipp_avg_decon_median=mean(decon_median,na.rm=TRUE),
                     Lhipp_avg_decon_sd=mean(decon_sd,na.rm=TRUE)) %>% ungroup()
  avg_Rhipp <- HCdata2 %>% group_by(id,run,trial,evt_time) %>% 
    dplyr::summarize(Rhipp_avg_decon_mean=mean(decon_mean,na.rm=TRUE),
                     Rhipp_avg_decon_median=mean(decon_median,na.rm=TRUE),
                     Rhipp_avg_decon_sd=mean(decon_sd,na.rm=TRUE)) %>% ungroup()
  #Rename the columns for merging
  avg_Lhipp <- avg_Lhipp %>% dplyr::rename(decon_mean=Lhipp_avg_decon_mean,
                                    decon_median=Lhipp_avg_decon_median,
                                    decon_sd=Lhipp_avg_decon_sd)
  avg_Rhipp <- avg_Rhipp %>% dplyr::rename(decon_mean=Rhipp_avg_decon_mean,
                                    decon_median=Rhipp_avg_decon_median,
                                    decon_sd=Rhipp_avg_decon_sd)
  #Add atlas values for HC data before combining
  avg_Lhipp$atlas_value <- 500
  avg_Rhipp$atlas_value <- 501
  #Merge HC data into subcortical_df 
  subcortical_df <- rbind(subcortical_df, avg_Lhipp)
  subcortical_df <- rbind(subcortical_df, avg_Rhipp)
  rm(HCdata1)
  rm(HCdata2)
  #Fix ID codings
  subcortical_df$id <- gsub(x = subcortical_df$id, pattern='881224',replacement='431224')
  subcortical_df$id <- gsub(x = subcortical_df$id, pattern='881230',replacement='431230')
  subcortical_df <- subcortical_df
  #fix how run is coded
  subcortical_df <- subcortical_df %>% mutate(run = case_when(
    str_detect(run, regex("run1", ignore_case=TRUE)) ~ 1,
    str_detect(run, regex("run2", ignore_case=TRUE)) ~ 2))
  #Add run_trial to medusa_df
  subcortical_df <- subcortical_df %>% mutate(trial = as.vector(trial)) %>% mutate(run_trial = case_when(
    trial < 121 ~ trial,
    trial > 120 ~ trial-120L)) # the L tells it that I want it to be an integer vector not a double vector
  #Filter data +/- 3 sec from event
  subcortical_df <- subcortical_df %>% filter(evt_time > -3 & evt_time < 3) 
  #Create labels for 12 slice HC data
  cobraHC_labels <- data.frame(roi_num=c(500,501),hemi=c('L','R'),network=c('Hippocampus_cobra','Hippocampus_cobra'),
                               subregion=c('Hippocampus_cobra','Hippocampus_cobra'),x=c(0,0),y=c(0,0),z=c(0,0),
                               MNI_Glasser_HCP_v1.0=c('Hippocampus_cobra','Hippocampus_cobra'),
                               Brainnetome_1.0=c('Hippocampus_cobra','Hippocampus_cobra'),
                               CA_ML_18_MNI=c('Hippocampus_cobra','Hippocampus_cobra'),
                               atlas_value=c(500,501)) 
  labels <- rbind(labels,cobraHC_labels)
  #Merge labels with medusa data
  subcortical_df <- subcortical_df %>% inner_join(labels, by="atlas_value") #add labels to subcortical_df 
  #Code the region groupings
  subcortical_df <- subcortical_df %>% mutate(group = as.factor(case_when(
    str_detect(subregion, regex("Anterior Putamen", ignore_case=TRUE)) ~ "Striatum",
    str_detect(subregion, regex("Caudate Tail and Lateral Putamen", ignore_case=TRUE)) ~ "Striatum",
    str_detect(subregion, regex("Caudate Head", ignore_case=TRUE)) ~ "Striatum",
    str_detect(subregion, regex("Ventral Striatum", ignore_case=TRUE)) ~ "Striatum",
    str_detect(subregion, regex("Posterior Putamen", ignore_case=TRUE)) ~ "Striatum",
    str_detect(subregion, regex("BLA", ignore_case=TRUE)) ~ "Amygdala",
    str_detect(subregion, regex("CMN", ignore_case=TRUE)) ~ "Amygdala",
    str_detect(subregion, regex("Hippocampus:", ignore_case=TRUE)) ~ "Hippocampus_schaeffer",
    #If want to code anterior vs posterior HC for the schaeffer data:
    #str_detect(subregion, regex("anterior", ignore_case=TRUE)) ~ "AH",
    #str_detect(subregion, regex("posterior", ignore_case=TRUE)) ~ "PH",
    str_detect(subregion, regex("Thalamus", ignore_case=TRUE)) ~ "Thalamus",
    str_detect(subregion, regex("Hippocampus_cobra", ignore_case=TRUE)) ~ "Hippocampus_cobra")))
  #Merge trial_df and medusa data
  Q <- merge(trial_df, subcortical_df, by = c("id", "run", "run_trial")) %>% dplyr::arrange("id","run","run_trial","evt_time")
  #Remove cerebellum data
  Q <- Q %>% filter(!subregion=="Cerebellum")
  Q$decon_mean[Q$evt_time > Q$iti_ideal] = NA
}
#Filter based on subject group you are running
if (subjs_to_run=="patients_only") {
  Q <- Q %>% filter(!Group == "Controls")
  Q$Group <- droplevels(Q$Group)
} else if (subjs_to_run=="controls_only") {
  Q <- Q %>% filter(Group == "Controls")
  Q$Group <- droplevels(Q$Group)
}
#Filter out networks, if desired
if (select_networks==TRUE) {
  Q <- Q %>% filter(network %in% network_subset)
  Q$network <- plyr::revalue(Q$network, c("Hippocampus_cobra" = "Hippocampus"))
}


  
########## Prepare final dataset that will be input into mixed_by models ########## 
#Censor out data that bleeds into adjacent trials
if(align_to == "response"){
  Q$decon_mean[Q$evt_time > Q$iti_ideal] = NA #censors next trial
  Q$decon_mean[Q$evt_time < -(Q$rt_csv + Q$iti_prev)] = NA # censors prior trial
} else if(align_to == "clock"){
  Q$decon_mean[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA # censors next trial
  Q$decon_mean[Q$evt_time < -Q$iti_prev] = NA # censors prior trial (edited) 
}
#Scale demographics variables
Q$female <- relevel(as.factor(Q$sex),ref="F")
Q$age <- scale(Q$age)
Q$wtar <- scale(Q$wtar)
Q$exit <- scale(Q$exit)
Q$education_yrs <- scale(Q$education_yrs)
Q$pid5_psychoticism <- scale(Q$pid5_psychoticism)
Q$uppsp_total <- scale(Q$uppsp_total)
Q$uppsp_negative_urgency <- scale(Q$uppsp_negative_urgency)
Q$uppsp_positive_urgency <- scale(Q$uppsp_positive_urgency)
Q$uppsp_lack_of_premeditation <- scale(Q$uppsp_lack_of_premeditation)
Q$uppsp_lack_of_perseveration <- scale(Q$uppsp_lack_of_perseveration)
Q$neo_neuroticism <- scale(Q$neo_neuroticism)
Q$drs_total <- scale(Q$drs_total)
Q$cirsg <- scale(Q$cirsg)
Q$ham_hamtotal_17items <- scale(Q$ham_hamtotal_17items)
Q$beckhopelessness <- scale(Q$beckhopelessness)
Q$athf <- scale(Q$athf)
Q$WSLS_rs <- scale(Q$WSLS_rs)
Q$LS_rs <- scale(Q$LS_rs)
Q <- Q %>% dplyr::rename(group_leth = groupLeth)
#Code groupLeth with controls as reference
Q$group_leth_c <- factor(Q$group_leth, level = c("Controls", "Depressed", "Ideators", "LL_Attempters", "HL_Attempters"))
Q$group_leth_LL <- factor(Q$group_leth, level = c("LL_Attempters",  "HL_Attempters", "Controls", "Depressed", "Ideators"))
Q$Group_d <- factor(Q$Group, level = c("Depressed","Controls","Ideators", "Attempters"))
Q$Group_i <-factor(Q$Group, level = c("Ideators", "Controls", "Depressed", "Attempters"))
#Fix duplicate trial entries; confirmed they are same with sum(Q$trial.x != Q$trial.y)
Q$trial <- Q$trial.x 
Q <- Q %>% select(-c(trial.x,trial.y))
#Filter out first 10 trials, if desired (to match what Andrew did)
if (filter_first_10_trials==TRUE) {
  Q <- Q %>% filter(trial>10)
}
#Added from Andrew (4/6/2024); run_trial0 goes from 1-40. Need this to run model that matches Andrew's models
Q <- Q %>% mutate(run_trial0 = case_when(trial <= 40 ~ trial, 
                                         trial > 40 & trial <= 80 ~ trial-40,
                                         trial > 80 & trial <=120 ~ trial-80, 
                                         trial > 120 & trial <=160 ~ trial-120,
                                         trial > 160 & trial <=200 ~ trial-160,
                                         trial > 200 & trial <=240 ~ trial-200))
Q <- Q %>% mutate(run_trial0_c = run_trial0-floor(run_trial0/40.5),
                  run_trial0_neg_inv = -(1 / run_trial0_c) * 100,
                  run_trial0_neg_inv_sc = as.vector(scale(run_trial0_neg_inv)))



########## Run mixed_by models ########## 
#Loop through all of the models
for (model_counter in 1:length(decode_formula)) {    
  print(model_counter)
  if (structure=="network") {
    splits = c('evt_time','network') #will do a separate MLM for each group in splits; evt_time is the TR; collapse across both hemispheres of subregions
    if(region_to_run == "subcortical") {
      dir <- paste0(subcortical_cache_dir,'/Network_CombinedHemis')
      dir.create(file.path(subcortical_cache_dir, '/Network_CombinedHemis'), showWarnings = FALSE)
    } else if(region_to_run == "cortical") {
      dir <- paste0(cortical_cache_dir,'/Network_CombinedHemis')
      dir.create(file.path(cortical_cache_dir, '/Network_CombinedHemis'), showWarnings = FALSE)
    }
    setwd(dir)
    df0 <- decode_formula[[model_counter]]
    print(df0)
    ddf <- mixed_by(Q, outcomes = "rt_csv", rhs_model_formulae = df0 , split_on = splits, return_models=TRUE,
                    padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3, #ncores put 1 less than total cores on your machine or else it will slow down
                    tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                    emmeans_spec = list(
                      RT = list(outcome='rt_csv', model_name='model1', 
                                specs=formula(~rt_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2),rt_lag_sc=c(-2,-1,0,1,2))),
                      Vmax = list(outcome='rt_csv', model_name='model1', 
                                  specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2),rt_vmax_lag_sc=c(-2,-1,0,1,2))),
                      RTxO = list(outcome='rt_csv',model_name='model1',
                                  specs=formula(~rt_lag_sc:last_outcome:subj_level_rand_slope), at=list(subj_level_rand_slope=c(-2,-1,0,1,2),rt_lag_sc=c(-2,-1,0,1,2)))        
                      
                    ),
                    emtrends_spec = list(
                      RT = list(outcome='rt_csv', model_name='model1', var='rt_lag_sc', 
                                specs=formula(~rt_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2))),
                      Vmax = list(outcome='rt_csv', model_name='model1', var='rt_vmax_lag_sc', 
                                  specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2))),
                      RTxO = list(outcome='rt_csv',model_name='model1',var='rt_lag_sc',
                                  specs=formula(~rt_lag_sc:last_outcome:subj_level_rand_slope), at=list(subj_level_rand_slope=c(-2,-1,0,1,2)))
                      
                    )
                    )
  }
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  if(align_to == "response"){
    if(region_to_run == "subcortical") {save(ddf,file=paste0(curr_date,'-subcortical-response-','model_',model_counter,output_file_specifier,'.Rdata')) # removed i (i is for if you are running several models)
    } else if (region_to_run == "cortical") {save(ddf,file=paste0(curr_date,'-vmpfc-response-','model_',model_counter,output_file_specifier,'.Rdata'))} # removed i (i is for if you are running several models)
  } else if(align_to == "clock"){
    if(region_to_run == "subcortical") {save(ddf,file=paste0(curr_date,'-subcortical-clock-','model_',model_counter,output_file_specifier,'.Rdata')) # removed i (i is for if you are running several models)
    } else if (region_to_run == "cortical") {save(ddf,file=paste0(curr_date,'-vmpfc-clock-','model_',model_counter,output_file_specifier,'.Rdata'))} # removed i (i is for if you are running several models)
  }
}



#Now run plot_mixed_by_AI.R; need to make clock_data_list.txt and response_data_list.txt before plotting