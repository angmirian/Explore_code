library(tidyverse)
library(psych)
library(MplusAutomation)
library(ggcorrplot)
library(dependlab)
library(ggplot2)
library(viridis)
library(cowplot)

rm(list=ls())

#this script can be used to apply network labels to extracted betas that were sig following FWE correction, and then generate corrplots to look at associations between brain regions in each network w/ antag and the random slopes of trustee/reciprocity

setwd("~/Documents/RESEARCH/ANALYSES/fMRI/April2023")
fmri_dir <- "~/Documents/RESEARCH/ANALYSES/fMRI/Dec2022/extracted_values/12Dec2022/"
bad_ids=c(207224, 210374, 210701, 213708, 216806, 220024, 221698, 440149, 219392, 216845, 440311, 440336) # AI - added all excludes
#labels <- read_csv("region_labels_244.csv") %>% mutate(roi_num = as.numeric(roi_num)) %>% inner_join(read_csv("region_lookup_244.csv"), by = "roi_num")
split_hipp_amyg = TRUE #if want to split amygdala and hippocampus for ROI analyses

align = "feedback" #option clock, feedback
region = "both" #options are vPFC, sub, and both
group = "pts" # options are pts and all (don't actually need to do anything for all)

#read in 17-network labels
map_717 <- read_csv("~/Documents/RESEARCH/ANALYSES/fMRI/Oct2022/FromTim_25Oct2022/Schaefer_200_7_17_mapping.csv") %>% 
  dplyr::select(-subregion7, -x, -y, -z, -hemi) #the subregion labels from the 244-parcel csv are more descriptive, so we'll use those since the two columns don't match (despite being conceptually the same). We don't need the coordinates, because those are provided in the pe_betas

# we need to add the extra 44 parcels. These are in the 7network labels, so we'll read those in and merge the extra 44 with our 17-network labels. 
labels <- read_csv(paste0("region_labels_244.csv")) %>% filter(!is.na(subregion)) %>% rename(roi_num7 = roi_num, network7 = network, subregion7 = subregion) #%>% dplyr::select(-hemi)
labels$roi_num7 <- as.numeric(labels$roi_num7)

#merge the 7 and 17 network labels
final_labels <- full_join(map_717, labels)

# Load subject data
sub_df <- readRDS("explore_n146.rds") %>% mutate(id = registration_redcapid)
pid5_df <- read.csv("~/Documents/RESEARCH/ANALYSES/Explore_behavioral/PID5/pid5_facets.csv")
sub_df <- left_join(sub_df,pid5_df,by="id")

#filter sub_df if just want pts
if (group == "pts") {
  sub_df <- sub_df %>% filter(!Group=="Controls")
} 

#demos <- sub_df %>%
#  dplyr::select(-id) %>% # ascending id
#  dplyr::rename(id=registration_redcapid,groupLeth=group_leth,sex=gender) %>% # real study id
#  group_by(id) %>%
#  dplyr::filter(row_number() == 1) %>%
#  ungroup() %>% 
#  #dplyr::select(id, age, education_yrs, group_leth, starts_with("pid5"), starts_with("ipip")) %>% 
#  dplyr::select(id, age, sex, education_yrs, wtar, exit, mmse, Group, Group_a, groupLeth, age_at_first_attempt, time_since_last_attempt, max_ideation, max_intent, names(pid5_df), starts_with("uppsp"), starts_with("neo"), starts_with("bis")) %>% 
#  droplevels() %>%
#  dplyr::mutate(id = as.integer(id), registration_edu = as.numeric(education_yrs))

#recode pid5 to start with "pid5"
sub_df <- sub_df %>%
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
sub_df<-sub_df %>% dplyr::mutate(pid5_negativeaffect = (pid5_anxiousness+ pid5_emotional_lability+ pid5_separation_insecurity)/3)
sub_df<-sub_df %>% dplyr::mutate(pid5_detachment = (pid5_withdrawal+ pid5_anhedonia+ pid5_intimacy_avoidance)/3)
sub_df<-sub_df %>% dplyr::mutate(pid5_antagonism = (pid5_manipulativeness+ pid5_deceitfulness+ pid5_grandiosity)/3)
sub_df<-sub_df %>% dplyr::mutate(pid5_disinhibition = (pid5_irresponsibility+ pid5_impulsivity+ pid5_distractibility)/3)
sub_df<-sub_df %>% dplyr::mutate(pid5_psychoticism = (pid5_unusual_beliefs_experiences+ pid5_eccentricity+ pid5_perceptual_dysregulation)/3)

# read in V_max betas
#setwd(file.path(fmri_dir, "L1m-v_max_wi")) #clock-aligned
#setwd(file.path(fmri_dir, "L1m-v_max_wi_lead")) #Jan 2023, response aligned

if (align=="clock") {
  setwd(file.path(fmri_dir, "L1m-v_max_wi")) #clock-aligned
  vmax <- read_csv("transformed_Schaefer_244_final_3.125mm_cope_l2.csv.gz") %>% filter(l2_model == "intercept" & l1_cope_name == "EV_v_max_wi" & !(id %in% bad_ids)) %>% rename(roi_num7 = mask_value) 
} else if (align=="feedback") {
  setwd(file.path(fmri_dir, "L1m-v_max_wi_lead")) #Jan 2023, feedback aligned
  vmax <- read_csv("transformed_Schaefer_244_final_3.125mm_cope_l2.csv.gz") %>% filter(l2_model == "intercept" & l1_cope_name == "EV_v_max_wi_lead" & !(id %in% bad_ids)) %>% rename(roi_num7 = mask_value) 
}

#merge labels with extracted betas
#vmax <- read_csv("transformed_Schaefer_244_final_3.125mm_cope_l2.csv.gz") %>% filter(l2_model == "intercept" & l1_cope_name == "EV_v_max" & !(id %in% bad_ids))
#vmax <- read_csv("transformed_Schaefer_244_final_3.125mm_cope_l2.csv.gz") %>% filter(l2_model == "intercept" & l1_cope_name == "EV_v_max_wi" & !(id %in% bad_ids)) %>% rename(roi_num7 = mask_value) 
#vmax <- read_csv("transformed_Schaefer_244_final_3.125mm_cope_l2.csv.gz") %>% filter(l2_model == "intercept" & l1_cope_name == "EV_v_max_wi_lead" & !(id %in% bad_ids)) %>% rename(roi_num7 = mask_value) 

#merge labels with extracted betas
#pe_betas <- left_join(pe_betas, final_labels, by = "roi_num7") %>% filter(l1_cope_name == "pe.overall") #I have confirmed the mask value is in roi7 units because 
#vmax <- vmax %>% rename(roi_num7 = mask_value) %>% select(-x, -y, -z) %>% filter(l1_cope_name == "pe.overall")
vmax <- vmax %>% left_join(final_labels, by = "roi_num7")

#vmax <- vmax %>% mutate(roi_num7 = mask_value) %>% select(-x, -y, -z) %>% inner_join(labels, by = "roi_num7") %>% inner_join(sub_df, by = "id")

# label missing networks as (amygdala, hippocampus, thalamus) - added option to split amygdala (2/26/2024)
if (split_hipp_amyg){
  labels$network[str_detect(labels$subregion, regex("Hippocampus", ignore_case=TRUE))] <- "hippocampus"
  labels$network[str_detect(labels$subregion, regex("BLA", ignore_case=TRUE))] <- "amygdala"
  labels$network[str_detect(labels$subregion, regex("CMN", ignore_case=TRUE))] <- "amygdala"
} else {labels$network[is.na(labels$network)] <- "amyg_hipp_thal"}

setwd("~/Documents/RESEARCH/ANALYSES/fMRI/April2023")

#filter to drop the visual network
#vmax <- vmax %>% filter(network7 != "Vis")

#this is all an attempt to get individual labels for each parcel 
vmax <- vmax %>% mutate(plot_label = sub("Focus point:\\s+", "", MNI_Glasser_HCP_v1.0, perl=TRUE)) %>% #default to glasser
  mutate(plot_label = sub("\\* Within \\d mm: ", "", plot_label, perl=TRUE)) %>% #remove extraneous characters 
  unite(subregion7, c(subregion7, hemi), sep = "_") %>% #to make subregion labels distinct, label them w/ hemisphere
  mutate(plot_label = if_else(is.na(plot_label), subregion7, plot_label)) %>% #use subregion labels for the extra 44 parcels
 mutate(plot_label = if_else(plot_label == "L_Area_dorsal_23_a+b", CA_ML_18_MNI, plot_label)) %>% #ca_ml_18 to resolve duplicates
  mutate(plot_label = if_else(plot_label == "R_Area_PGs", CA_ML_18_MNI, plot_label))

#remove residual language introduced from mni labels
vmax <- vmax %>% mutate(plot_label = sub("Focus point:\\s+", "", plot_label, perl=TRUE)) %>% 
  mutate(plot_label = sub("\\* Within \\d mm: ", "", plot_label, perl=TRUE))

#Filter to just vPFC and subcortical ROIs
if (region=="vPFC"){
  vPFC <- vmax %>% 
    filter(roi_num7 %in% c(84,86,88,89,192,194,65,66,67,170,171,55,56,159,160,161,191))
  vmax <- vPFC
} else if (region =="sub") {
  sub <- vmax %>% 
    filter(roi_num7 %in% c(201,202,203,204,223,224,225,226,227,228,229,230,215,220,216,221,217,222))
  vmax <- sub
} else if (region == "both") {
  both <- vmax %>%
    filter(roi_num7 %in% c(84,86,88,89,192,194,65,66,67,170,171,55,56,159,160,161,191,201,202,203,204,223,224,225,226,227,228,229,230,215,220,216,221,217,222)) 
  vmax <- both
}

#lots of duplicate labels still, so append the roi_number
vmax <- vmax %>% unite(plot_label, c(plot_label, roi_num7))

#wide form 
vmax.wide <- vmax %>% pivot_wider(id_cols = id, names_from = plot_label, values_from = value) 

## Not sure what's going on from here down - messaged Tim to get an example of what's in these files
#loop over the two models
#mname <- c("Documents/GitHub/bsocial_trust_fmri/trust_mplus/wm_rs.out", "Documents/GitHub/bsocial_trust_fmri/trust_mplus/wm_rs_1trustee.out", "Documents/GitHub/bsocial_trust_fmri/trust_mplus/wm_rs_int.out")
#slopevars <- c("REX.Mean", "RGOOD.Mean", "RBAD.Mean", "RRECP.Mean", "RRECP_EX.Mean")
#filename <- c("pe7network_traits_designvars_me_corrplot.pdf", "pe7network_traits_designvars_1t_corrplot.pdf", "pe7network_traits_designvars_interaction_corrplot.pdf")

#for (i in 1:3) {
  #read in model
#  wm_rs <- readModels(paste0(os, mname[i]))
  
  ## save the savedata from the mplus output to an object
#  slopes <- wm_rs$savedata %>% 
#    dplyr::select(ID, any_of(slopevars)) %>% 
#    group_by(ID) %>% slice_head() %>% ungroup() %>%
#    rename(id = ID)
  
vmax.rs <- left_join(vmax.wide, sub_df, by = "id")
  
#filename <- c("pe7network_traits_designvars_COGNITION_corrplot.pdf")

  #corrplots
#  fscores <- read_csv(paste0(os, "Documents/GitHub/bsocial_trust_fmri/agree_factors_6-21-22.csv"))
  #for plots, it will help to have unique variance for callousness / narcissism
#  m1 <- lm(callous2f ~ narc2f, data = fscores)
#  fscores$call_res <- m1$residuals
#  m2 <- lm(narc2f ~ callous2f, data = fscores)
#  fscores$narc_res <- m2$residuals
  
#  pe.rs <- left_join(pe.rs, fscores, by = "id")
  
################## COMMENTED CODE FROM TIM

#this part may or may not be necessary for you. I wanted to generate 7 heatmaps -- one for each functional network. pe_betas is a long-form dataframe where I had the parcel labels and which networks they belonged to. Here, I'm just getting the network assignments and then merging them into the object I created above (targetcors)
#note that in targetcors, the withvar variable will be = to the names of your parcels
map2roi17 <- vmax %>% dplyr::select(network7, plot_label)
map2roi17 <- unique(map2roi17) %>% rename(withvar = plot_label)

#targetcors <- left_join(targetcors, map2roi17)

#write.csv(targetcors, "~/Documents/GitHub/bsocial_trust_fmri/tc_corrmatrix.csv", row.names = F)

#this is a loop to generate the heatmap plots
#networks <- c("Cont", "Default", "DorsAttn", "Limbic", "SalVentAttn", "SomMot")
#out <- list()
#for (q in networks) {
#  net_df <- targetcors %>% filter(network7 == paste0(q))
#  p <- ggplot(net_df, aes(targetvar, withvar, fill= r)) + ggtitle(paste0(q)) +
#    geom_tile() + geom_text(aes(label = r)) + scale_fill_viridis(discrete=FALSE, limits=c(-.3, .3)) + theme_cowplot() + theme(axis.text.x=element_text(angle=90,hjust=1))
#  out[[q]] <- p
#}

#net_plots <- cowplot::plot_grid(plotlist = out)

#cowplot::save_plot(paste0(gd, "Hallquist_Dombrovski/Social_Trust/plots/", filename[i]), net_plots, base_height = 15, base_width = 20)


########################


#source the function
source(paste0('~/Documents/RESEARCH/ANALYSES/fMRI/scripts/From_Tim_6April2023/taacor.R'), echo=TRUE)

#create vector of column names for behavioral variables
#slope_tar <- names(slopes[,2:ncol(slopes)])"
#slope_tar <- c("exit", "wtar", "uppsp_negative_urgency","uppsp_lack_of_premeditation", "uppsp_lack_of_perseveration","uppsp_positive_urgency","neo_neuroticism","neo_extraversion","neo_openness","neo_aggreableness","neo_conscientiousness","pid5_negativeaffect","pid5_detachment","pid5_antagonism","pid5_disinhibition","pid5_psychoticism")

#vmax.rs <- vmax.rs %>% dplyr::select(id:Cerebellum_R_244 | exit | wtar | uppsp_negative_urgency:neo_conscientiousness | pid5_negativeaffect:pid5_psychoticism)
vmax.rs <- vmax.rs %>% dplyr::select(id:'Hippocampus: anterior head_R_230' | exit | wtar | uppsp_negative_urgency:neo_conscientiousness | pid5_negativeaffect:pid5_psychoticism)
vmax.rs <- vmax.rs %>% mutate_if(is.character,as.numeric)



#pe.rs is a wide-form dataframe with my betas, my behavioral variables, and my personality variables
# the function will make a limited heatmap. The target argument specifies which variables should be on the x-axis, everything else will be on the y-axis
# absrmin tells it to suppress any correlations < absolute value of .15
# targetcors <- taacor(pe.rs %>% dplyr::select(-id), target = c("antag1f", "narc_res", "call_res", slope_tar), omit = c(), digits = 2, absrmin = .15)

#filename <- c("vmax_clock_7network_traits_designvars_COGNITION_corrplot.pdf")
#filename <- c("vmax_clock_7network_traits_designvars_UPPSP_corrplot.pdf")
#filename <- c("vmax_clock_7network_traits_designvars_NEO_corrplot.pdf")
#filename <- c("vmax_clock_7network_traits_designvars_PID5_corrplot.pdf")
#filename <- c("vmax_response_7network_traits_designvars_COGNITION_corrplot.pdf")
#filename <- c("vmax_response_7network_traits_designvars_UPPSP_corrplot.pdf")
#filename <- c("vmax_response_7network_traits_designvars_NEO_corrplot.pdf")
filename <- c("vmax_response_7network_traits_designvars_PID5_corrplot.pdf")

#slope_tar <- c("exit", "wtar", "uppsp_negative_urgency","uppsp_lack_of_premeditation", "uppsp_lack_of_perseveration","uppsp_positive_urgency","neo_neuroticism","neo_extraversion","neo_openness","neo_aggreableness","neo_conscientiousness","pid5_negativeaffect","pid5_detachment","pid5_antagonism","pid5_disinhibition","pid5_psychoticism")
#targetcors <- taacor(vmax.rs %>% dplyr::select(-id), target = c("exit", "wtar"), omit = c(), digits = 2, absrmin = .15)
#targetcors <- taacor(vmax.rs %>% dplyr::select(-id), target = c("uppsp_negative_urgency","uppsp_lack_of_premeditation", "uppsp_lack_of_perseveration","uppsp_positive_urgency"), omit = c(), digits = 2, absrmin = .15)
#targetcors <- taacor(vmax.rs %>% dplyr::select(-id), target = c("neo_neuroticism","neo_extraversion","neo_openness","neo_aggreableness","neo_conscientiousness"), omit = c(), digits = 2, absrmin = .15)
targetcors <- taacor(vmax.rs %>% dplyr::select(-id), target = c("pid5_negativeaffect","pid5_detachment","pid5_antagonism","pid5_disinhibition","pid5_psychoticism"), omit = c(), digits = 2, absrmin = .15)

  
  #map2roi17 <- pe_betas %>% dplyr::select(network7, plot_label)
  map2roi17 <- vmax %>% dplyr::select(network7, plot_label)
  map2roi17 <- unique(map2roi17) %>% rename(withvar = plot_label)
  
  targetcors <- left_join(targetcors, map2roi17)
  
  #write.csv(targetcors, "~/Documents/GitHub/bsocial_trust_fmri/tc_corrmatrix.csv", row.names = F)
  #networks <- c("Cont", "Default", "DorsAttn", "Limbic", "SalVentAttn", "SomMot","amyg_hipp_thal")
  networks <- c("Cont", "Default", "Limbic", "amyg_hipp_thal")
  out <- list()
  for (q in networks) {
    net_df <- targetcors %>% filter(network7 == paste0(q))
    p <- ggplot(net_df, aes(targetvar, withvar, fill= r)) + ggtitle(paste0(q)) +
      geom_tile() + geom_text(aes(label = r), size=16) + scale_fill_viridis(discrete=FALSE, limits=c(-.3, .3)) + theme_cowplot() + 
      theme(#axis.text = element_text(size=30),
            legend.title = element_text(size=40, face = "bold"), #change legend title font size
            legend.text = element_text(size=20),  
           axis.text.x = element_text(angle=90,hjust=1, size=40),
            axis.text.y = element_text(size=40),
            plot.title = element_text(size = 50, face = "bold"),
            axis.title = element_text(size=45))
    out[[q]] <- p
  }
  net_plots <- cowplot::plot_grid(plotlist = out)
  #cowplot::save_plot(paste0("~/Documents/RESEARCH/ANALYSES/fMRI/April2023/plots/", filename[i]), net_plots, base_height = 15, base_width = 20)
  cowplot::save_plot(paste0("~/Documents/RESEARCH/ANALYSES/fMRI/April2023/plots/", filename), net_plots, base_height = 30, base_width = 40)
  
  #cowplot::save_plot(paste0(gd, "Hallquist_Dombrovski/Social_Trust/plots/", filename[i]), net_plots, base_height = 15, base_width = 20)
#}
