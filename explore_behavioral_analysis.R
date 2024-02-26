# Script to try to replicate Aliona's behavioral findings in Explore sample
#Aliona's results for personality measures: https://dnpl.slite.com/app/docs/ykGR_9PMw9sNsO
#Angela Ianni
#Feb 14, 2023
#July 12, 2023 - added code to get demographic information for Table 1 of paper
#Feb 14, 2024 - ungrouped before scaling condition_trial_inv

#rm(list=ls())

library(tidyverse)
library(lme4)
library(emmeans)
library(parallel)
library(broom)
library(broom.mixed)
library(RColorBrewer)

setwd("~/Documents/RESEARCH/ANALYSES/Explore_behavioral")

exclude_controls=0 # set to 1 to exclude controls from analysis

# behavioral analyses
#source("~/Documents/RESEARCH/ANALYSES/Explore_behavioral/scripts/get_trial_data.R") #from https://github.com/UNCDEPENdLab/clock_analysis/blob/master/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R
#trial_df <- get_trial_data("~/Documents/RESEARCH/ANALYSES/temporal_instrumental_agent/clock_task/vba_fmri", dataset = "explore") 
#saveRDS(trial_df, file = "trial_df_02142023.rds")
#trial_df <- readRDS("trial_df_02142023.rds")
trial_df <- readRDS("trial_df_03282023.rds")
trial_df$id <- as.integer(sub("_1", "", trial_df$id))

#Load subject's data
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
    #dplyr::select(id, age, education_yrs, group_leth, starts_with("pid5"), starts_with("ipip")) %>% 
    dplyr::select(id, age, sex, race, education_yrs, wtar, exit, mmse, beckhopelessness, ham_hamtotal_17items, cirsg, drs_total, Group, Group_a, 
                  groupLeth, age_at_first_attempt, time_since_last_attempt, max_ideation, max_intent, 
                  athf, opioid, antipsychotic, sedhyp, 
                  anxiety_dx, suds_dx,
                  names(pid5_df), starts_with("uppsp"), starts_with("neo"), starts_with("bis")) %>% 
    droplevels() %>%
    dplyr::mutate(id = as.integer(id), registration_edu = as.numeric(education_yrs))

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

#Remove sensitivity analyses - behavior flatliners
#demos_orig <- demos
#demos <- demos_orig
excludes = c(221292, 221637, 440263) # For Aliona's paper
demos <- demos %>% filter(!id %in% excludes)

#Remove excludes - imaging QC (n=12)
#excludes = c(207224, 210374, 210701, 213708, 216806, 216845, 219392, 220024, 221698, 440149, 440311, 440336)
#demos <- demos %>% filter(!id %in% excludes)

#Remove all excludes and sensitivity analyses subjects - imaging QC
#excludes = c(207224, 210374, 210701, 213708, 216806, 216845, 219392, 220024, 221698, 440149, 440311, 440336, 211253, 431100, 440243, 440254, 440263, 440392, 221292, 221637)
#demos <- demos %>% filter(!id %in% excludes)

#trial_df <- trial_df %>% filter(!id %in% excludes1)
trial_df <- trial_df %>% filter(!id %in% excludes)


#Get demographic information
#need to add race, dementia rating scale, physical illness burden, hamilton rating scale for depression score, 
#beck hopelessness scale score, antidepressant exposure, lifetime substance use (no), lifetime anxiety (no)
summary(demos$Group)
summ <- demos %>% group_by(Group) %>%
    summarise(group_n = n(),
              males_n = sum(sex=="M"),
              males_percent = 100*(males_n/group_n),
              age_mean = mean(age),
              age_sd = sd(age),
              edu_yrs_mean = mean(education_yrs),
              edu_yrs_sd = sd(education_yrs), 
              premorbid_iq_n = sum(!is.na(wtar)),
              premorbid_iq = mean(wtar, na.rm=TRUE),
              premorbid_iq_sd = sd(wtar, na.rm=TRUE),
              exit_n = sum(!is.na(exit)),
              exit_mean = mean(exit, na.rm=TRUE),
              exit_sd = sd(exit, na.rm=TRUE),
              drs_mean = mean(drs_total, na.rm=TRUE),
              drs_sd = sd(drs_total, na.rm=TRUE),
              cirsg_mean = mean(cirsg, na.rm=TRUE),
              cirsg_sd = sd(cirsg, na.rm=TRUE),
              ham_mean = mean(ham_hamtotal_17items, na.rm=TRUE),
              ham_sd = sd(ham_hamtotal_17items, na.rm=TRUE),
              beckhopelessness_mean = mean(beckhopelessness, na.rm=TRUE),
              beckhopelessness_sd = sd(beckhopelessness, na.rm=TRUE))
#summ has the data that is contained in Table 1
summ_all_subjs <- demos %>%
    summarise(group_n = n(),
              males_n = sum(sex=="M"),
              males_percent = 100*(males_n/group_n),
              age_mean = mean(age),
              age_sd = sd(age),
              age_min = min(age),
              age_max = max(age),
              edu_yrs_mean = mean(education_yrs),
              edu_yrs_sd = sd(education_yrs), 
              premorbid_iq_n = sum(!is.na(wtar)),
              premorbid_iq = mean(wtar, na.rm=TRUE),
              premorbid_iq_sd = sd(wtar, na.rm=TRUE),
              exit_n = sum(!is.na(exit)),
              exit_mean = mean(exit, na.rm=TRUE),
              exit_sd = sd(exit, na.rm=TRUE),
              drs_mean = mean(drs_total, na.rm=TRUE),
              drs_sd = sd(drs_total, na.rm=TRUE),
              cirsg_mean = mean(cirsg, na.rm=TRUE),
              cirsg_sd = sd(cirsg, na.rm=TRUE),
              ham_mean = mean(ham_hamtotal_17items, na.rm=TRUE),
              ham_sd = sd(ham_hamtotal_17items, na.rm=TRUE),
              beckhopelessness_mean = mean(beckhopelessness, na.rm=TRUE),
              beckhopelessness_sd = sd(beckhopelessness, na.rm=TRUE))

#Group by controls vs. all patients
demos <- demos %>% mutate(simple_groups = case_when(
    str_detect(Group, regex("Controls", ignore_case=TRUE)) ~ "Controls",
    str_detect(Group, regex("Depressed", ignore_case=TRUE)) ~ "Patients",
    str_detect(Group, regex("Ideators", ignore_case=TRUE)) ~ "Patients",
    str_detect(Group, regex("Attempters", ignore_case=TRUE)) ~ "Patients",
))
summary(demos$simple_groups)
demos %>% group_by(simple_groups) %>%
    summarise(group_n = n(),
              males_n = sum(sex=="M"),
              males_percent = 100*(males_n/group_n),
              age_mean = mean(age),
              age_sd = sd(age),
              edu_yrs_mean = mean(education_yrs),
              edu_yrs_sd = sd(education_yrs), 
              premorbid_iq_n = sum(!is.na(wtar)),
              premorbid_iq = mean(wtar, na.rm=TRUE),
              premorbid_iq_sd = sd(wtar, na.rm=TRUE),
              exit_n = sum(!is.na(exit)),
              exit_mean = mean(exit, na.rm=TRUE),
              exit_sd = sd(exit, na.rm=TRUE),
              drs_mean = mean(drs_total, na.rm=TRUE),
              drs_sd = sd(drs_total, na.rm=TRUE),
              cirsg_mean = mean(cirsg, na.rm=TRUE),
              cirsg_sd = sd(cirsg, na.rm=TRUE),
              ham_mean = mean(ham_hamtotal_17items, na.rm=TRUE),
              ham_sd = sd(ham_hamtotal_17items, na.rm=TRUE),
              beckhopelessness_mean = mean(beckhopelessness, na.rm=TRUE),
              beckhopelessness_sd = sd(beckhopelessness, na.rm=TRUE)) %>% ungroup()


#Run statistical tests - for simple groups (all pts vs. controls - ACNP 2023 abstract)
sex.df = table(demos$sex, demos$simple_groups)
print(sex.df)
print(chisq.test(sex.df))
age_aov <- aov(age ~ simple_groups, data = demos)
summary(age_aov)
edu_aov <- aov(education_yrs ~ simple_groups, data = demos)
summary(edu_aov)
wtar_aov <- aov(wtar ~ simple_groups, data = demos)
summary(wtar_aov)
exit_aov <- aov(exit ~ simple_groups, data = demos)
summary(exit_aov)
cor.test(demos$age, demos$exit, use="pairwise.complete.obs")

race.df = table(demos$race, demos$simple_groups)
print(race.df)
print(chisq.test(race.df))

#Run statistical tests - for 4 groups
sex.df = table(demos$sex, demos$Group)
print(sex.df)
print(chisq.test(sex.df))
age_aov <- aov(age ~ Group, data = demos)
summary(age_aov)
edu_aov <- aov(education_yrs ~ Group, data = demos)
summary(edu_aov)
pairwise.t.test(demos$education_yrs, demos$Group, p.adj = "none") #post-hoc test
wtar_aov <- aov(wtar ~ Group, data = demos)
summary(wtar_aov)
pairwise.t.test(demos$wtar, demos$Group, p.adj = "none") #post-hoc test
exit_aov <- aov(exit ~ Group, data = demos)
summary(exit_aov)
race.df = table(demos$race, demos$Group)
print(race.df)
print(chisq.test(race.df))

drs_aov <- aov(drs_total ~ Group, data = demos)
summary(drs_aov)
pairwise.t.test(demos$drs_total, demos$Group, p.adj = "none") #post-hoc test
cirsg_aov <- aov(cirsg ~ Group, data = demos) #physcial illness burden
summary(cirsg_aov)
pairwise.t.test(demos$cirsg, demos$Group, p.adj = "none") #post-hoc test
ham_aov <- aov(ham_hamtotal_17items ~ Group, data = demos)
summary(ham_aov)
pairwise.t.test(demos$ham_hamtotal_17items, demos$Group, p.adj = "none") #post-hoc test
beckhopelessness_aov <- aov(beckhopelessness ~ Group, data = demos)
summary(beckhopelessness_aov)
pairwise.t.test(demos$beckhopelessness, demos$Group, p.adj = "none") #post-hoc test


pt_demos <- demos %>% filter(Group != "Controls")
antidep.df = table(pt_demos$athf, pt_demos$Group)
print(antidep.df)
antidep_aov <- aov(athf ~ Group, data = pt_demos)
summary(antidep_aov)
sudx.df = table(pt_demos$suds_dx, pt_demos$Group)
print(sudx.df)
print(chisq.test(pt_demos$suds_dx, pt_demos$Group))
anx.df = table(pt_demos$anxiety_dx, pt_demos$Group)
print(anx.df)
print(chisq.test(pt_demos$anxiety_dx, pt_demos$Group))
pairwise.chisq.test(pt_demos$anxiety_dx, pt_demos$Group, p.adjust.method = "none")

#summary(demos$groupLeth)
#summary(demos$age)
#sd(demos$age,na.rm=TRUE)
#table(demos$sex)
#summary(demos$education_yrs)
#sd(demos$education_yrs,na.rm=TRUE)
#summary(demos$wtar)
#sd(demos$wtar,na.rm=TRUE)
#summary(demos$mmse)
#sd(demos$mmse,na.rm=TRUE)
#summary(demos$age_at_first_attempt)
#sd(demos$age_at_first_attempt,na.rm=TRUE)
#summary(demos$time_since_last_attempt)
#sd(demos$time_since_last_attempt,na.rm=TRUE)
#summary(demos$max_ideation)
#sd(demos$max_ideation,na.rm=TRUE)
#summary(demos$max_intent)
#sd(demos$max_intent,na.rm=TRUE)

trial_df <- trial_df %>%
    left_join(demos, by="id") %>%
    mutate(reward_lag_rec = if_else(reward_lag=="omission", -0.5, 0.5),
           reward_lag2 = !omission_lag2,
           reward_lag3 = !omission_lag3,
           reward_lag4 = !omission_lag4)

#Code the condition trial (more similar to bsocial run_trial)
#trial_df <- trial_df %>% 
#    group_by(id,run_number) %>% 
#    arrange(id, run_number, trial) %>% 
#    mutate(condition_trial = run_trial-floor(run_trial/40.5)*40,                                                                                          
#           condition_trial_neg_inv = -1000 / condition_trial,
#           condition_trial_neg_inv_sc = as.vector(scale(condition_trial_neg_inv)), 
#           ev_sc = scale(ev),
#           ev_lag_sc = scale(lag(ev))) %>% ungroup()

#Should ungroup before scaling condition_trial_neg_inv (Per Alex on 2/14/2024)
trial_df <- trial_df %>% 
    group_by(id,run_number) %>% 
    arrange(id, run_number, trial) %>% 
    mutate(condition_trial = run_trial-floor(run_trial/40.5)*40,                                                                                          
           condition_trial_neg_inv = -1000 / condition_trial) %>% ungroup
trial_df <- trial_df %>% mutate(condition_trial_neg_inv_sc = as.vector(scale(condition_trial_neg_inv)), 
           ev_sc = scale(ev),
           ev_lag_sc = scale(lag(ev)))


#If want to exclude controls
if (exclude_controls==1) {
    trial_df <- trial_df[trial_df$Group!="Controls",]
    trial_df$Group <- droplevels(trial_df$Group)
}






###May 22, 2023 - plots for Aliona's paper
##Figure 2
trial_df$groupLeth <- factor(trial_df$groupLeth, level = c("HL_Attempters", "LL_Attempters", "Ideators", "Depressed", "Controls"))
#Get info on how many pts per group
#trial_df[trial_df$groupLeth=="HL_Attempters",]
dfx3_filt <- trial_df

#Plots for Aliona - Updated 9/29/23
##plotting rt_csv in seconds by trial; panel by contingency; AI - change run trial to condition trial
b1<-ggplot(dfx3_filt, aes(condition_trial, rt_csv, color = rewFunc)) + geom_smooth(method = "gam", linewidth = 3) +
    labs(x="Trial", y="Response time (seconds)", color = "Contingency") +
    scale_color_manual(values = c("#8B4513", "black"))+
                                   theme(axis.title=element_text(size=20)) +
                                       theme(axis.text = element_text(size = 16)) +
                                       theme(legend.title = element_text(size = 20)) +
                                       theme(legend.text = element_text(size = 16))
print(b1)
#With df control: formula = y~splines::ns(x,2))
b2<-ggplot(dfx3_filt, aes(condition_trial, rt_swing, color = rewFunc)) +  geom_smooth(method="gam", formula = y~splines::ns(x,5), linewidth = 3) +
   labs(x="Trial", y="Change in RT (seconds)", color = "Contingency") +
   scale_color_manual(values = c("#8B4513", "black"))+
                                  theme(axis.title=element_text(size=20)) +
                                      theme(axis.text = element_text(size = 16)) +
                                      theme(legend.title = element_text(size = 20)) +
                                      theme(legend.text = element_text(size = 16))
print(b2)
#7/27/23 - look at AGE EFFECTS on E/E behavior
dfx3_filt_wide <- dfx3_filt  %>% distinct_at(vars(id), .keep_all=TRUE)
###values for plotting approximated based on descriptives and visual inspection of distributions
psych::describe(dfx3_filt_wide$age) #range 49-80, mean 62.01 (median 61), sd is 6.8 (mean-2*sd = 48.4 -> 49, mean+2*sd = 75.6 -> 76; mean-1*sd = 55, mean+1*sd = 69)
#11/28/2023 - added stuff for ACNP 2023
#replace rt_vmax with ev (n.s.)
#m100 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(age) + ev_lag_sc*condition_trial_neg_inv_sc*scale(age) + rt_lag_sc*condition_trial_neg_inv_sc*scale(age)+ (1|id/run), dfx3_filt)

#Original model
m100 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(age) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(age) + rt_lag_sc*condition_trial_neg_inv_sc*scale(age)+ (1|id/run), dfx3_filt)
car::Anova(m100,'3')
summary(m100)
####PLOTTING  rt lag*reward lag*Trait  at approximately +-2 Gmds (2 SD is 49, 76; 1 SD is 55, 69)
em100 <- as_tibble(emmeans::emtrends(m100, var = "rt_lag_sc", specs = c("reward_lag_rec", "age"), at=list(age=c(49,76)), options = list() ))
em100$reward_lag_rec <- factor(em100$reward_lag_rec, levels = c("-0.5", "0.5"), labels = c("Reward omission", "Reward"))
# Convert age to factor
em100$age <- factor(em100$age, levels = c("49", "76"), labels = c("49y", "76y"))
brewer_palette_tmp <- brewer.pal(n = 3, name = "Set2")
#Plot RT swings
ggplot(em100, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=age)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +
    ylab("RT swings\n Small <---------> Large") +
    theme(
        axis.title = element_text(size = 14),  # Increase axis title font size
        axis.text = element_text(size = 12, color = "grey10"),  # Increase axis tick label font size and set color
        panel.grid.major.x = element_blank(),  # Remove major x-axis grid lines
        legend.text = element_text(size = 12),  # Increase legend label font size
        legend.title = element_text(size = 14)  # Increase legend title font size
    ) +
    scale_x_discrete(name = "") +
    scale_y_reverse() + labs (color = "Age") +
    scale_color_brewer(palette = "Set2")
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em100b <- as_tibble(emmeans::emtrends(m100, var = "rt_vmax_lag_sc", specs = c("age", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), age = c(49,76)), options = list() ))
em100b$age <- factor(em100b$age, levels = c("49", "76"), labels = c("49y", "76y"))
ggplot(em100b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=age)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(color="Age") +
    #labs(color = "Lethality") +
    theme(axis.text=element_text(size=8.5, color="grey10")) + 
    scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40")) + 
    scale_color_brewer(palette = "Set2")

#7/27/23 - look at EXIT EFFECTS on E/E behavior
#dfx3_filt_wide <- dfx3_filt  %>% distinct_at(vars(id), .keep_all=TRUE)
###values for plotting approximated based on descriptives and visual inspection of distributions
psych::describe(dfx3_filt_wide$exit) #range 0-16, mean 5.81 (median 5), sd is 2.96 (mean-2*sd = 0, mean+2*sd = 12)
m101 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(exit) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(exit) + rt_lag_sc*condition_trial_neg_inv_sc*scale(exit)+ (1|id/run), dfx3_filt)
car::Anova(m101,'3')
summary(m101)
####PLOTTING  rt lag*reward lag*Trait  at approximately +-2 Gmds
em101 <- as_tibble(emmeans::emtrends(m101, var = "rt_lag_sc", specs = c("reward_lag_rec", "exit"), at=list(exit=c(0,12)), options = list() ))
em101$reward_lag_rec <- factor(em101$reward_lag_rec, levels = c("-0.5", "0.5"), labels = c("Reward omission", "Reward"))
# Convert exit to factor
em101$exit <- factor(em101$exit, levels = c("0", "12"), labels = c("Low", "High"))
brewer_palette_tmp <- brewer.pal(n = 3, name = "Set2")
#Plot RT swings
ggplot(em101, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=exit)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +
    ylab("RT swings\n Small <---------> Large") +
    theme(
        axis.title = element_text(size = 14),  # Increase axis title font size
        axis.text = element_text(size = 12, color = "grey10"),  # Increase axis tick label font size and set color
        panel.grid.major.x = element_blank(),  # Remove major x-axis grid lines
        legend.text = element_text(size = 12),  # Increase legend label font size
        legend.title = element_text(size = 14)  # Increase legend title font size
    ) +
    scale_x_discrete(name = "") +
    scale_y_reverse() + labs (color = "Exit") +
    scale_color_brewer(palette = "Set2")
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em101b <- as_tibble(emmeans::emtrends(m101, var = "rt_vmax_lag_sc", specs = c("exit", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), exit = c(0,12)), options = list() ))
ggplot(em101b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=exit)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))

#Test out Working Memory Strategy Utilization
m102 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(age) + rt_lag2_sc*reward_lag2*scale(age) + rt_lag3_sc*reward_lag3*scale(age) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(age) + rt_lag_sc*condition_trial_neg_inv_sc*scale(age)+ (1|id/run), dfx3_filt)
car::Anova(m102,'3')
summary(m102)
m103 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(exit) + rt_lag2_sc*reward_lag2*scale(exit) + rt_lag3_sc*reward_lag3*scale(exit) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(exit) + rt_lag_sc*condition_trial_neg_inv_sc*scale(exit)+ (1|id/run), dfx3_filt)
car::Anova(m103,'3')
summary(m103)

## For Aliona's paper
excludes = c(221292, 221637, 440263) # behavioral flatliners
demos <- demos %>% filter(!id %in% excludes)
trial_df <- trial_df %>% filter(!id %in% excludes)
trial_df$groupLeth <- factor(trial_df$groupLeth, level = c("HL_Attempters", "LL_Attempters", "Ideators", "Depressed", "Controls"))
dfx3_filt <- trial_df
m3a <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*groupLeth + rt_vmax_lag_sc*condition_trial_neg_inv*groupLeth + rt_lag_sc*condition_trial_neg_inv*groupLeth + (1|id/run), dfx3_filt)
car::Anova(m3a,'3')
summary(m3a)

#pair_m3a.emm <- emmeans(m3a, pairwise ~ reward_lag_rec*groupLeth)
#pairs(pair_m3a.emm, adjust="fdr")

#pair_em2n <- emmeans(m2n, pairwise~Group_a | network)
#pairs(pair_em2n, adjust="fdr")

####PLOTTING RT SWINGS
em3a <- as_tibble(emmeans::emtrends(m3a, var = "rt_lag_sc", specs = c("reward_lag_rec", "groupLeth")))
em3a$reward_lag_rec <- factor(em3a$reward_lag_rec, levels = c("-0.5", "0.5"), labels = c("Reward omission", "Reward"))
#specify the number of colors (n) and the palette type (type) to generate the desired color scheme
#Here the brewer.pal() function generates a palette with 4 colors from the "Set1" palette.
library(RColorBrewer)
# Define the ColorBrewer palette
brewer_palette_tmp <- brewer.pal(n = 5, name = "Set1")
brewer_palette <- brewer_palette_tmp
brewer_palette[3] <- brewer_palette_tmp[5]
brewer_palette[4] <- brewer_palette_tmp[3]
brewer_palette[5] <- brewer_palette_tmp[4]
#brewer_palette <- brewer.pal(n = 4, name = "Set1")
# Update the color values in scale_color_manual()
ggplot(em3a, aes(x = reward_lag_rec, y = rt_lag_sc.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = groupLeth)) +
    geom_point(position = position_dodge(width = .6), size = 2.5) + 
    geom_errorbar(position = position_dodge(width = 0.6), width = 0.4, size = 0.9) + 
    theme_bw(base_size = 12) +
    ylab("RT swings\n Small <---------> Large") +
    theme(
        axis.title = element_text(size = 14),  # Increase axis title font size
        axis.text = element_text(size = 12, color = "grey10"),  # Increase axis tick label font size and set color
        panel.grid.major.x = element_blank(),  # Remove major x-axis grid lines
        legend.text = element_text(size = 12),  # Increase legend label font size
        legend.title = element_text(size = 14)  # Increase legend title font size
    ) +
    scale_x_discrete(name = "") +
    scale_y_reverse() +
    scale_color_manual(name = "Group",
                       values = brewer_palette,
                       labels = c("High-Lethality", "Low-Lethality", "Ideators", "Depressed", "Controls"))
#Personality:
##filter to get patients only dataset
brewer_palette_tmp <- brewer.pal(n = 3, name = "Set2")
dfx3_filt_subset_pts <- dfx3_filt
dfx3_filt_subset_pts <- dfx3_filt_subset_pts[trial_df$Group!="Controls",]
dfx3_filt_subset_pts$Group <- droplevels(dfx3_filt_subset_pts$Group)
#dfx3_filt_subset_BPD <- dfx3_filt %>% filter(groupLeth!="HC") 
dfx3_filt_subset_pts_wide <- dfx3_filt_subset_pts  %>% distinct_at(vars(id), .keep_all=TRUE)
###values for plotting approximated based on descriptives and visual inspection of distributions
#describe(dfx3_filt_subset_BPD_wide$ipip_a) #range 35-115, Gmd 16.5 (32.74 is 2 SDs); mean 85.72;
psych::describe(dfx3_filt_subset_pts_wide$neo_aggreableness) #range 33-55, mean 45.7 (median 46), sd is 4.85 (mean-2*sd = 36, mean+2*sd = 55.4)
m13 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(neo_aggreableness) + rt_vmax_lag_sc*condition_trial_neg_inv*scale(neo_aggreableness) + rt_lag_sc*condition_trial_neg_inv*scale(neo_aggreableness)+ (1|id/run), dfx3_filt_subset_pts)
####PLOTTING  rt lag*reward lag*Trait  at approximately +-2 Gmds
emm13a <- as_tibble(emmeans::emtrends(m13, var = "rt_lag_sc", specs = c("reward_lag_rec", "neo_aggreableness"), at=list(neo_aggreableness=c(36,55)), options = list() ))
emm13a$reward_lag_rec <- factor(emm13a$reward_lag_rec, levels = c("-0.5", "0.5"), labels = c("Reward omission", "Reward"))
# Convert neo_aggreableness to factor
emm13a$neo_aggreableness <- factor(emm13a$neo_aggreableness, levels = c("36", "55"), labels = c("Low", "High"))
#emm13a <- emm13a %>% rename(neo_agreeableness=neo_aggreableness)
ggplot(emm13a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=neo_aggreableness)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +
    ylab("RT swings\n Small <---------> Large") +
    theme(
        axis.title = element_text(size = 14),  # Increase axis title font size
        axis.text = element_text(size = 12, color = "grey10"),  # Increase axis tick label font size and set color
        panel.grid.major.x = element_blank(),  # Remove major x-axis grid lines
        legend.text = element_text(size = 12),  # Increase legend label font size
        legend.title = element_text(size = 14)  # Increase legend title font size
    ) +
    scale_x_discrete(name = "") +
    scale_y_reverse() + labs (color = "Agreeableness") +
    scale_color_brewer(palette = "Set2")
#negative affect
m10 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_negativeaffect) + rt_vmax_lag_sc*condition_trial_neg_inv*scale(pid5_negativeaffect) + rt_lag_sc*condition_trial_neg_inv*scale(pid5_negativeaffect)+ (1|id/run), dfx3_filt_subset_pts)
#describe(dfx3_filt_subset_BPD_wide$pid5_negativeaffect) #range 0.583-2.83, Gmd 0.6935 (1.387 is 2 SDs); mean 2.053;
psych::describe(dfx3_filt_subset_pts_wide$pid5_negativeaffect) #range 0-2.42, mean = 1.05 (median 1.11), sd 0.62, mean-2*sd=-0.19 (will use 0), mean+2*sd=2.29
emm10a <- as_tibble(emmeans::emtrends(m10, var = "rt_lag_sc", specs = c("reward_lag_rec", "pid5_negativeaffect"), at=list(pid5_negativeaffect=c(0,2.3)), options = list() ))
emm10a$reward_lag_rec <- factor(emm10a$reward_lag_rec, levels = c("-0.5", "0.5"), labels = c("Reward omission", "Reward"))
# Convert ipip_a to factor
emm10a$pid5_negativeaffect <- factor(emm10a$pid5_negativeaffect, levels = c("0", "2.3"), labels = c("Low", "High"))
####PLOTTING  rt lag*reward lag*Trait
ggplot(emm10a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_negativeaffect)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +
    ylab("RT swings\n Small <---------> Large") +
    theme(
        axis.title = element_text(size = 14),  # Increase axis title font size
        axis.text = element_text(size = 12, color = "grey10"),  # Increase axis tick label font size and set color
        panel.grid.major.x = element_blank(),  # Remove major x-axis grid lines
        legend.text = element_text(size = 12),  # Increase legend label font size
        legend.title = element_text(size = 14)  # Increase legend title font size
    ) +
    scale_x_discrete(name = "") +
    scale_y_reverse() + labs (color = "Negative affect") +
    scale_color_brewer(palette = "Set2")


###Compare behavior across groups
m1 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*Group + rt_vmax_lag_sc*condition_trial_neg_inv_sc*Group +
                rt_lag_sc*condition_trial_neg_inv_sc*Group + (1|id/run), trial_df)
car::Anova(m1,'3')
em1 <- as_tibble(emmeans::emtrends(m1, var = "rt_lag_sc", specs = c("reward_lag_rec", "Group") ))
ggplot(em1, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=Group)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) +  
    scale_y_reverse()
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em1b <- as_tibble(emmeans::emtrends(m1, var = "rt_vmax_lag_sc", specs = c("Group", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial)), options = list() ))
ggplot(em1b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=Group)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))


###Compare lethality groups
trial_df$groupLeth <- factor(trial_df$groupLeth, level = c("Controls", "Depressed", "Ideators", "LL_Attempters", "HL_Attempters"))
#original model - this is what is in Aliona's paper (I think - as of 9/11/2023)
m2 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*groupLeth + rt_vmax_lag_sc*condition_trial_neg_inv_sc*groupLeth +
                rt_lag_sc*condition_trial_neg_inv_sc*groupLeth+ (1|id/run), trial_df)
#add random slopes
m2 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*groupLeth + rt_vmax_lag_sc*condition_trial_neg_inv_sc*groupLeth +
               rt_lag_sc*condition_trial_neg_inv_sc*groupLeth + (rt_lag_sc + reward_lag_rec|id) + (1|id/run), trial_df)
#add interaction random slopes
m2 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*groupLeth + rt_vmax_lag_sc*condition_trial_neg_inv_sc*groupLeth +
               rt_lag_sc*condition_trial_neg_inv_sc*groupLeth + (rt_lag_sc*reward_lag_rec|id) + (1|id/run), trial_df)
#Code from Aliona (4/11/23)
##Rt lag as random
m2 <- lmer(scale(rt_csv) ~ rt_lag_sc*reward_lag_rec*groupLeth + rt_vmax_lag_sc*condition_trial_neg_inv_sc*groupLeth +
               rt_lag_sc*condition_trial_neg_inv_sc*groupLeth + (rt_lag_sc|id) + (1|id/run), trial_df,
           control = lmerControl(optimizer ="Nelder_Mead")) # as from https://stats.stackexchange.com/questions/242109/model-failed-to-converge-warning-in-lmer
while (any(grepl("failed to converge", m2@optinfo$conv$lme4$messages) )) {
    print(m2@optinfo$conv$lme4$conv)
    ss <- getME(m2,c("theta","fixef"))
    m2 <- update(m2, start=ss)}
#while (any(grepl("failed to converge", m2@optinfo$conv$lme4$messages) )) {
#    print(m2@optinfo$conv$lme4$conv)
#    ss <- getME(m2,c("theta","fixef"))
#    m2 <- update(m2, start=ss)}
##reward lag rec as random - model failed to converge: degenerate Hessian with 1 negative eigenvalues
m2 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*groupLeth + rt_vmax_lag_sc*condition_trial_neg_inv_sc*groupLeth +
               rt_lag_sc*condition_trial_neg_inv_sc*groupLeth + (reward_lag_rec|id) + (1|id/run), trial_df,
           control = lmerControl(optimizer ="Nelder_Mead")) # as from https://stats.stackexchange.com/questions/242109/model-failed-to-converge-warning-in-lmer
while (any(grepl("failed to converge", m2@optinfo$conv$lme4$messages) )) {
    print(m2@optinfo$conv$lme4$conv)
    ss <- getME(m2,c("theta","fixef"))
    m2 <- update(m2, start=ss)}
#while (any(grepl("failed to converge", m2@optinfo$conv$lme4$messages) )) {
#    print(m2@optinfo$conv$lme4$conv)
#    ss <- getME(m2,c("theta","fixef"))
#    m2 <- update(m2, start=ss)}

car::Anova(m2, type=3)
#summary(m2)
em2 <- as_tibble(emmeans::emtrends(m2, var = "rt_lag_sc", specs = c("reward_lag_rec", "groupLeth") ))
p <- ggplot(em2, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=groupLeth)) +
    #scale_x_discrete(name ="Outcome of Previous Trial", labels=c("Omission", "Reward")) +
    geom_point(position = position_dodge(width = .7), size=4) + 
    geom_errorbar(position = position_dodge(width=0.7), width=0.8, size=1.1) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    xlab("Outcome of Previous Trial") + 
    theme(axis.text.y=element_blank(), panel.grid.major.y=element_blank(),
         axis.text.x=element_blank(), panel.grid.major.x=element_blank(),
          #axis.text.x=element_text(face="bold",size=14, color="grey10"),
          axis.title=element_text(face="bold",size=20),
          legend.text=element_text(size=14),
          legend.title=element_text(face="bold",size=20)) +
          #axis.text.x = element_text(face="bold",  size=14)) +
    scale_y_reverse() +
    scale_color_discrete(name = "Group", labels = c("Controls", "Depressed", "Ideators", "Attempters - Low Lethality", "Attempters - High Lethality")) #+
print(p)

# breaks=c("-0.5","0.5"),
#Plot exploitation - not significant
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em2b <- as_tibble(emmeans::emtrends(m2, var = "rt_vmax_lag_sc", specs = c("groupLeth", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial)), options = list() ))
ggplot(em2b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=groupLeth)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))


###Split by age of first attempt
m3 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(age_at_first_attempt) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(age_at_first_attempt) +
               rt_lag_sc*condition_trial_neg_inv_sc*scale(age_at_first_attempt) + (1|id/run), trial_df)
car::Anova(m3, type=3)
summary(m3)
em3a <- as_tibble(emmeans::emtrends(m3, var = "rt_lag_sc", specs = c("reward_lag_rec", "age_at_first_attempt"), at = list(age_at_first_attempt = c(20, 60))))
ggplot(em3a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=age_at_first_attempt)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) +  
    scale_y_reverse()
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em3b <- as_tibble(emmeans::emtrends(m3, var = "rt_vmax_lag_sc", specs = c("age_at_first_attempt", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), age_at_first_attempt = c(20, 60)), options = list() ))
ggplot(em3b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=age_at_first_attempt)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))

###Time since last attempt - not significant
m4 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(time_since_last_attempt) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(time_since_last_attempt) +
               rt_lag_sc*condition_trial_neg_inv_sc*scale(time_since_last_attempt) + (1|id/run), trial_df)
car::Anova(m4, type=3)

###Max intent
m5 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(max_intent) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(max_intent) +
               rt_lag_sc*condition_trial_neg_inv_sc*scale(max_intent) + (1|id/run), trial_df)
car::Anova(m5, type=3)
summary(m5)
avg_trait <- mean(trial_df$max_intent, na.rm=TRUE)
sd_trait <- sd(trial_df$max_intent, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em5a <- as_tibble(emmeans::emtrends(m5, var = "rt_lag_sc", specs = c("reward_lag_rec", "max_intent"), at = list(max_intent = c(low_trait, high_trait))))
ggplot(em5a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=max_intent)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) +  
    scale_y_reverse()
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em5b <- as_tibble(emmeans::emtrends(m5, var = "rt_vmax_lag_sc", specs = c("max_intent", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), max_intent = c(low_trait, high_trait)), options = list() ))
ggplot(em5b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=max_intent)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))







### Compare behavior between groups 
library(tidyverse)
library(ggpubr)
library(rstatix)
#summary_stats_score <- trial_df %>% group_by(id, Group) %>% get_summary_stats(score_csv, type = "mean_sd")
summary_stats_score <- trial_df %>% group_by(id, simple_groups) %>% get_summary_stats(score_csv, type = "mean_sd")
summary_stats_score_individual <- trial_df %>% group_by(id) %>% get_summary_stats(score_csv, type = "mean_sd")
summary_stats_score_individual_IEV <- trial_df %>% group_by(id) %>% filter(rewFunc=="IEV") %>% get_summary_stats(score_csv, type = "mean_sd")
summary_stats_score_individual_DEV <- trial_df %>% group_by(id) %>% filter(rewFunc=="DEV") %>% get_summary_stats(score_csv, type = "mean_sd")
score_individual_demos.df <- left_join(demos, summary_stats_score_individual, by="id")
score_individual_IEV_demos.df <- left_join(demos, summary_stats_score_individual_IEV, by="id")
score_individual_DEV_demos.df <- left_join(demos, summary_stats_score_individual_DEV, by="id")

#Check model for IEV and DEV separately
trial_df_IEV <-trial_df %>% filter(rewFunc=="IEV")
trial_df_DEV <- trial_df %>% filter(rewFunc=="DEV")
trial_df_early <- trial_df %>% filter(condition_trial<20)
trial_df_late <- trial_df %>% filter(condition_trial>=20)

m100_all <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(age) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(age) + rt_lag_sc*condition_trial_neg_inv_sc*scale(age)+ (1|id/run), trial_df)
m100_IEV <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(age) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(age) + rt_lag_sc*condition_trial_neg_inv_sc*scale(age)+ (1|id/run), trial_df_IEV)
m100_DEV <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(age) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(age) + rt_lag_sc*condition_trial_neg_inv_sc*scale(age)+ (1|id/run), trial_df_DEV)
m100_early <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(age) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(age) + rt_lag_sc*condition_trial_neg_inv_sc*scale(age)+ (1|id/run), trial_df_early)
m100_late <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(age) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(age) + rt_lag_sc*condition_trial_neg_inv_sc*scale(age)+ (1|id/run), trial_df_late)
car::Anova(m100_all, type=3)
summary(m100_all)
car::Anova(m100_IEV, type=3)
car::Anova(m100_DEV, type=3)
summary(m100_IEV)
summary(m100_DEV)
car::Anova(m100_early, type=3)
summary(m100_early)
car::Anova(m100_late, type=3)
summary(m100_late)

#Plots for ALL trials/conditions
#plot at 1st and 3rd quartiles (from summary(demos$age))
em100_all <- as_tibble(emmeans::emtrends(m100_all, var = "rt_lag_sc", specs = c("reward_lag_rec", "age"), at=list(age=c(49,76)), options = list() )) #1 SD = 57, 67; mean=62.03, SD=6.82; median = 61, min 49, max 80; 2SD = 48.3 -> 49, 75.6 -> 76
em100_all$age <- factor(em100_all$age, levels = c("49", "76"), labels = c("mean-2*sd (49y)", "mean+2*sd (76y)"))
brewer_palette_tmp <- brewer.pal(n = 3, name = "Set2")
#Plot RT swings
ggplot(em100_all, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=age)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +
    ylab("RT swings\n Small <---------> Large") +
    theme(
        axis.title = element_text(size = 14),  # Increase axis title font size
        axis.text = element_text(size = 12, color = "grey10"),  # Increase axis tick label font size and set color
        panel.grid.major.x = element_blank(),  # Remove major x-axis grid lines
        legend.text = element_text(size = 12),  # Increase legend label font size
        legend.title = element_text(size = 14)  # Increase legend title font size
    ) +
    scale_x_discrete(name = "") +
    scale_y_reverse() + labs (color = "Age") +
    scale_color_brewer(palette = "Set2")
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em100b_all <- as_tibble(emmeans::emtrends(m100_all, var = "rt_vmax_lag_sc", specs = c("age", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), age = c(49,76)), options = list() ))
em100b_all$age <- factor(em100b_all$age, levels = c("49", "76"), labels = c("mean-2*sd (49y)", "mean+2*sd (76y)"))
ggplot(em100b_all, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=age)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(color="Age") +
    #labs(color = "Lethality") +
    theme(axis.text=element_text(size=8.5, color="grey10")) + 
    scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40")) + 
    scale_color_brewer(palette = "Set2")
#Plot exploration by trial effect
em100c_all <- as_tibble(emmeans::emtrends(m100_all, var = "rt_lag_sc", specs = c("age", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), age = c(49,76)), options = list() ))
em100c_all$age <- factor(em100c_all$age, levels = c("49", "76"), labels = c("mean-2*sd (49y)", "mean+2*sd (76y)"))
ggplot(em100c_all, aes(x=condition_trial_neg_inv_sc, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=age)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("RT swings\n Small <---------> Large")  +
    labs(color="Age") +
    #labs(color = "Lethality") +
    theme(axis.text=element_text(size=8.5, color="grey10")) + 
    scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40")) + 
    scale_color_brewer(palette = "Set2")

###Plots for IEV ONLY
em100_IEV <- as_tibble(emmeans::emtrends(m100_IEV, var = "rt_lag_sc", specs = c("reward_lag_rec", "age"), at=list(age=c(49,76)), options = list() ))
em100_IEV$age <- factor(em100_IEV$age, levels = c("49", "76"), labels = c("mean-2*sd (49y)", "mean+2*sd (76y)"))
brewer_palette_tmp <- brewer.pal(n = 3, name = "Set2")
#Plot RT swings
ggplot(em100_IEV, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=age)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +
    ylab("RT swings\n Small <---------> Large") +
    theme(
        axis.title = element_text(size = 14),  # Increase axis title font size
        axis.text = element_text(size = 12, color = "grey10"),  # Increase axis tick label font size and set color
        panel.grid.major.x = element_blank(),  # Remove major x-axis grid lines
        legend.text = element_text(size = 12),  # Increase legend label font size
        legend.title = element_text(size = 14)  # Increase legend title font size
    ) +
    scale_x_discrete(name = "") +
    scale_y_reverse() + labs (color = "Age") +
    scale_color_brewer(palette = "Set2")
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em100b_IEV <- as_tibble(emmeans::emtrends(m100_IEV, var = "rt_vmax_lag_sc", specs = c("age", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), age = c(49,76)), options = list() ))
em100b_IEV$age <- factor(em100b_IEV$age, levels = c("49", "76"), labels = c("mean-2*sd (49y)", "mean+2*sd (76y)"))
ggplot(em100b_IEV, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=age)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(color="Age") +
    #labs(color = "Lethality") +
    theme(axis.text=element_text(size=8.5, color="grey10")) + 
    scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40")) + 
    scale_color_brewer(palette = "Set2")
#Plot exploration by trial effect
em100c_IEV <- as_tibble(emmeans::emtrends(m100_IEV, var = "rt_lag_sc", specs = c("age", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), age = c(49,76)), options = list() ))
em100c_IEV$age <- factor(em100c_IEV$age, levels = c("49", "76"), labels = c("mean-2*sd (49y)", "mean+2*sd (76y)"))
ggplot(em100c_IEV, aes(x=condition_trial_neg_inv_sc, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=age)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("RT swings\n Small <---------> Large")  +
    labs(color="Age") +
    #labs(color = "Lethality") +
    theme(axis.text=element_text(size=8.5, color="grey10")) + 
    scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40")) + 
    scale_color_brewer(palette = "Set2")


###Plots for DEV ONLY
em100_DEV <- as_tibble(emmeans::emtrends(m100_DEV, var = "rt_lag_sc", specs = c("reward_lag_rec", "age"), at=list(age=c(49,76)), options = list() ))
em100_DEV$age <- factor(em100_DEV$age, levels = c("49", "76"), labels = c("mean-2*sd (49y)", "mean+2*sd (76y)"))
brewer_palette_tmp <- brewer.pal(n = 3, name = "Set2")
#Plot RT swings
ggplot(em100_DEV, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=age)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +
    ylab("RT swings\n Small <---------> Large") +
    theme(
        axis.title = element_text(size = 14),  # Increase axis title font size
        axis.text = element_text(size = 12, color = "grey10"),  # Increase axis tick label font size and set color
        panel.grid.major.x = element_blank(),  # Remove major x-axis grid lines
        legend.text = element_text(size = 12),  # Increase legend label font size
        legend.title = element_text(size = 14)  # Increase legend title font size
    ) +
    scale_x_discrete(name = "") +
    scale_y_reverse() + labs (color = "Age") +
    scale_color_brewer(palette = "Set2")
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em100b_DEV <- as_tibble(emmeans::emtrends(m100_DEV, var = "rt_vmax_lag_sc", specs = c("age", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), age = c(49,76)), options = list() ))
em100b_DEV$age <- factor(em100b_DEV$age, levels = c("49", "76"), labels = c("mean-2*sd (49y)", "mean+2*sd (76y)"))
ggplot(em100b_DEV, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=age)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(color="Age") +
    #labs(color = "Lethality") +
    theme(axis.text=element_text(size=8.5, color="grey10")) + 
    scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40")) + 
    scale_color_brewer(palette = "Set2")
#Plot exploration by trial effect
em100c_DEV <- as_tibble(emmeans::emtrends(m100_DEV, var = "rt_lag_sc", specs = c("age", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), age = c(49,76)), options = list() ))
em100c_DEV$age <- factor(em100c_DEV$age, levels = c("49", "76"), labels = c("mean-2*sd (49y)", "mean+2*sd (76y)"))
ggplot(em100c_DEV, aes(x=condition_trial_neg_inv_sc, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=age)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("RT swings\n Small <---------> Large")  +
    labs(color="Age") +
    #labs(color = "Lethality") +
    theme(axis.text=element_text(size=8.5, color="grey10")) + 
    scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40")) + 
    scale_color_brewer(palette = "Set2")



#check for correlations
pear_whole <- cor.test(score_individual_demos.df$age, score_individual_demos.df$mean, use="pairwise.complete.obs")
spear_whole <- cor.test(score_individual_demos.df$age, score_individual_demos.df$mean, method="spearman", use="pairwise.complete.obs", exact = FALSE)
pear_IEV <- cor.test(score_individual_IEV_demos.df$age, score_individual_IEV_demos.df$mean, use="pairwise.complete.obs")
spear_IEV <- cor.test(score_individual_IEV_demos.df$age, score_individual_IEV_demos.df$mean, method="spearman", use="pairwise.complete.obs", exact = FALSE)
pear_DEV <- cor.test(score_individual_DEV_demos.df$age, score_individual_DEV_demos.df$mean, use="pairwise.complete.obs")
spear_DEV <- cor.test(score_individual_DEV_demos.df$age, score_individual_DEV_demos.df$mean, method="spearman", use="pairwise.complete.obs", exact = FALSE)

#check if groups differ in mean points earned
mean_score_aov <- aov(mean ~ simple_groups, data = score_individual_demos.df)
age_aov <- aov(age ~ simple_groups, data = score_individual_demos.df)



#combine into single dataframe
score_individual_demos.df <- score_individual_demos.df %>% rename(mean_score_whole=mean, sd_whole=sd)
score_individual_IEV_demos.df <- score_individual_IEV_demos.df %>% rename(mean_score_IEV=mean, sd_IEV=sd)
score_individual_DEV_demos.df <- score_individual_DEV_demos.df %>% rename(mean_score_DEV=mean, sd_DEV=sd)
score_individual_demos.df <- left_join(score_individual_demos.df,score_individual_IEV_demos.df,by="id")
score_individual_demos.df <- left_join(score_individual_demos.df,score_individual_DEV_demos.df,by="id")



#Plot age by mean score for different conditions
ggplot(data=score_individual_demos.df, aes(x=age)) + 
    geom_point(aes(y=mean_score_whole), color="black",size=3) + geom_smooth(aes(y=mean_score_whole), method=lm, color="black", size=3) +
    geom_point(aes(y=mean_score_IEV), color="blue",size=3) + geom_smooth(aes(y=mean_score_IEV), method=lm, color="blue", size=3) +
    geom_point(aes(y=mean_score_DEV), color="darkorange",size=3) + geom_smooth(aes(y=mean_score_DEV), method=lm, color="darkorange", size=3) +
    ylab("Mean Score per Trial") +
    xlab("Age (years)") +
    scale_fill_discrete(breaks=c("mean_score_whole","mean_score_IEV","mean_score_DEV")) +
    theme(legend.title = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          #panel.background = element_rect(fill = background_color),
          axis.title.y = element_text(margin=margin(r=6)),
          axis.title.x = element_text(margin=margin(t=6)))

ggplot(data=trial_df, aes(x=v_max, y=rt_vmax, color=ev)) + geom_point() +
    scale_color_gradient2(midpoint = 30, low = "blue", mid = "white", high = "red")
ggplot(data=trial_df, aes(x=rt_csv, y=rt_vmax, color=ev)) + geom_point() + geom_smooth(method=lm, color="yellow") +
    scale_color_gradient2(midpoint = 30, low = "blue", mid = "white", high = "red")
trial_df_RTltfive <- trial_df %>% filter(rt_csv<4.999)
ggplot(data=trial_df_RTltfive, aes(x=rt_csv, y=rt_vmax, color=ev)) + geom_point() + geom_smooth(method=lm, color="yellow") +
    scale_color_gradient2(midpoint = 30, low = "blue", mid = "white", high = "red")
cor.test(trial_df_RTltfive$rt_vmax, trial_df_RTltfive$rt_csv)
summary_stats_corr_rt_rtvmax <- trial_df_RTltfive %>% group_by(id) %>% 
    mutate(rt_rtvmax_corr = cor(rt_csv,rt_vmax,use="complete.obs"))
#summary_stats_corr_rt_rtvmax_demos.df <- left_join(demos,summary_stats_corr_rt_rtvmax,by="id") 
summary_stats_corr_rt_rtvmax_demos.df <- summary_stats_corr_rt_rtvmax %>% select(id,age,exit,rt_rtvmax_corr) %>% unique()
cor.test(summary_stats_corr_rt_rtvmax_demos.df$age,summary_stats_corr_rt_rtvmax_demos.df$rt_rtvmax_corr)





score_df <- trial_df %>% group_by(id) %>% mutate(total_score = sum(score_csv)) %>% select(id, Group, simple_groups, total_score, uppsp_total, uppsp_negative_urgency, uppsp_positive_urgency, uppsp_lack_of_perseveration, uppsp_lack_of_premeditation) %>% unique() %>% ungroup()
#bxp <- ggboxplot(summary_stats_score, x="Group", y="mean", color="Group", palatte = "jco")
bxp <- ggboxplot(summary_stats_score, x="simple_groups", y="mean", color="simple_groups", palatte = "jco")
#summary_stats_score <- trial_df %>% group_by(id, simple_groups) %>% get_summary_stats(score_csv, type = "mean_sd")
#score_df <- trial_df %>% group_by(id) %>% mutate(total_score = sum(score_csv)) %>% select(id, age, simple_groups, total_score, uppsp_total, uppsp_negative_urgency, uppsp_positive_urgency, uppsp_lack_of_perseveration, uppsp_lack_of_premeditation) %>% unique() %>% ungroup()
#bxp <- ggboxplot(summary_stats_score, x="simple_groups", y="mean", color="simple_groups", palatte = "jco")
bxp
#trial_df %>% group_by(id, Group) %>% identify_outliers(score_csv) #several - 219 x 192
#trial_df %>% group_by(id, Group) %>% shapiro_test(score_csv) #146x5
#ggqqplot(trial_df, "score_csv", ggtheme = theme_bw()) + facet_grid(trial~Group)
res.aov <- rstatix::anova_test(
    data = score_df, dv = total_score, wid = id,
    between = Groups)
#res.aov <- rstatix::anova_test(
#    data = score_df, dv = total_score, wid = id,
#    between = simple_groups)
get_anova_table(res.aov)

#10/30/2023 - check for impulsivity behavioral effects
cor.test(score_df$total_score, score_df$uppsp_total, use="pairwise.complete.obs")

#check for impulsivity subcategories
cor.test(score_df$total_score, score_df$uppsp_negative_urgency, use="pairwise.complete.obs")
cor.test(score_df$total_score, score_df$uppsp_positive_urgency, use="pairwise.complete.obs")
cor.test(score_df$total_score, score_df$uppsp_lack_of_perseveration, use="pairwise.complete.obs")
cor.test(score_df$total_score, score_df$uppsp_lack_of_premeditation, use="pairwise.complete.obs")

#Look for correlations between uppsp subscales
cor.test(score_df$uppsp_negative_urgency, score_df$uppsp_positive_urgency, use="pairwise.complete.obs")
cor.test(score_df$uppsp_negative_urgency, score_df$uppsp_lack_of_perseveration, use="pairwise.complete.obs")
cor.test(score_df$uppsp_negative_urgency, score_df$uppsp_lack_of_premeditation, use="pairwise.complete.obs")
cor.test(score_df$uppsp_positive_urgency, score_df$uppsp_lack_of_perseveration, use="pairwise.complete.obs")
cor.test(score_df$uppsp_positive_urgency, score_df$uppsp_lack_of_premeditation, use="pairwise.complete.obs")
cor.test(score_df$uppsp_lack_of_perseveration, score_df$uppsp_lack_of_premeditation, use="pairwise.complete.obs")

mydata <- score_df %>% select(uppsp_total, uppsp_negative_urgency, uppsp_positive_urgency, uppsp_lack_of_perseveration, uppsp_lack_of_premeditation)
cormat <- round(cor(mydata, use="pairwise.complete.obs"),2)
#cormat <- rcorr(as.matrix(mydata, use="pairwise.complete.obs"),type="pearson")
library(reshape2)
melted_cormat <- melt(cormat)
library(ggplot2)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    geom_text(aes(label = value)) +
    scale_fill_gradient2(low = "white", high = "red", 
                          limit = c(0,1), space = "Lab", 
                          name="Pearson\nCorrelation") +

    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))+
    coord_fixed()






#Do same for two groups (all pts vs. controls)
summary_stats_score_2groups <- trial_df %>% group_by(id, simple_groups, age) %>% get_summary_stats(score_csv, type = "mean_sd")
total_score.df <- trial_df %>% group_by(id) %>% select(id, age, simple_groups, total_score) %>% unique() %>% ungroup()
res.aov_2groups <- anova_test(
    data = trial_df, dv = score_csv, wid = id,
    between = simple_groups)
get_anova_table(res.aov_2groups)
#cor.test(summary_stats_score_2groups$age, summary_stats_score_2groups$mean)
cor.test(total_score.df$age, total_score.df$total_score)
res.aov_2groups_totalscore <- anova_test(
    data = total_score.df, dv = total_score, wid = id,
    between = simple_groups)
get_anova_table(res.aov_2groups_totalscore)
#Plot age effects on mean score by group
ggplot(data=summary_stats_score_2groups, aes(x=age, y=mean, color=simple_groups)) + geom_point() + geom_smooth(method=lm)
ggplot(data=total_score.df, aes(x=age, y=total_score, color=simple_groups)) + geom_point() + geom_smooth(method=lm)



#anova_group_score <- aov(score_csv~Group, data=trial_df_id)
#summary(anova_group_score) 
#tukey_group <- TukeyHSD(anova_group_score)
#plot(tukey_group,las=1)
#ggplot(trial_df_id, aes(x=Group, y=score_csv, group=Group)) + geom_point(cex=1.5, pch = 1.0, position = position_jitter(w=0.1, h=0))
#anova_groupleth_score <- aov(score_csv~groupLeth, data=trial_df_id)
#summary(anova_groupleth_score) 
#TukeyHSD(anova_groupleth_score)
res1.aov <- anova_test(
    data = trial_df, dv = score_csv, wid = id,
    between = groupLeth)
get_anova_table(res1.aov)
summary_stats_score1 <- trial_df %>% group_by(id, groupLeth) %>% get_summary_stats(score_csv, type = "mean_sd")
ggboxplot(summary_stats_score1, x="groupLeth", y="mean", color="groupLeth", palatte = "jco")

### Compare HL vs. LL attempters
anova_groupLeth_age <- aov(age~group_leth, data=sub_df)
summary(anova_groupLeth_age) #n.s.
anova_groupLeth_agefirstatt <- aov(age_at_first_attempt~group_leth, data=sub_df)
summary(anova_groupLeth_agefirstatt)
anova_groupLeth_maxintent <- aov(max_intent~group_leth, data=sub_df)
summary(anova_groupLeth_maxintent)
anova_groupLeth_wtar <- aov(wtar~group_leth, data=sub_df)
summary(anova_groupLeth_wtar)
TukeyHSD(anova_groupLeth_wtar)
anova_groupLeth_exit <- aov(exit~group_leth, data=sub_df)
summary(anova_groupLeth_exit)
anova_groupLeth_neoC <- aov(neo_conscientiousness~group_leth, data=sub_df)
summary(anova_groupLeth_neoC)
TukeyHSD(anova_groupLeth_neoC)
anova_groupLeth_neoN <- aov(neo_neuroticism~group_leth, data=sub_df)
summary(anova_groupLeth_neoN)
TukeyHSD(anova_groupLeth_neoN)
anova_groupLeth_neoE <- aov(neo_extraversion~group_leth, data=sub_df)
summary(anova_groupLeth_neoE)
TukeyHSD(anova_groupLeth_neoE)
anova_groupLeth_neoA <- aov(neo_aggreableness~group_leth, data=sub_df)
summary(anova_groupLeth_neoA)
TukeyHSD(anova_groupLeth_neoA)
anova_groupLeth_neoO <- aov(neo_openness~group_leth, data=sub_df)
summary(anova_groupLeth_neoO)
TukeyHSD(anova_groupLeth_neoO)
#UPPSP
anova_uppsp_total <- aov(uppsp_total~Group, data=demos)
summary(anova_uppsp_total)
TukeyHSD(anova_uppsp_total)
anova_groupLeth_uppsp_total <- aov(uppsp_total~groupLeth, data=demos)
summary(anova_groupLeth_uppsp_total)
TukeyHSD(anova_groupLeth_uppsp_total)

anova_groupLeth_uppspNU <- aov(uppsp_negative_urgency~group_leth, data=sub_df)
summary(anova_groupLeth_uppspNU)
TukeyHSD(anova_groupLeth_uppspNU)
anova_groupLeth_uppspLpremed <- aov(uppsp_lack_of_premeditation~group_leth, data=sub_df)
summary(anova_groupLeth_uppspLpremed)
TukeyHSD(anova_groupLeth_uppspLpremed)
anova_groupLeth_uppspLpersev <- aov(uppsp_lack_of_perseveration~group_leth, data=sub_df)
summary(anova_groupLeth_uppspLpersev)
TukeyHSD(anova_groupLeth_uppspLpersev)
anova_groupLeth_uppspPU <- aov(uppsp_positive_urgency~group_leth, data=sub_df)
summary(anova_groupLeth_uppspPU)
TukeyHSD(anova_groupLeth_uppspPU)
#neo
trial_df %>% dplyr::group_by(groupLeth) %>% dplyr::summarize(mean=mean(neo_conscientiousness, na.rm=TRUE))
trial_df %>% dplyr::group_by(groupLeth) %>% dplyr::summarize(mean=mean(neo_neuroticism, na.rm=TRUE))
trial_df %>% dplyr::group_by(groupLeth) %>% dplyr::summarize(mean=mean(neo_extraversion, na.rm=TRUE))
trial_df %>% dplyr::group_by(groupLeth) %>% dplyr::summarize(mean=mean(neo_aggreableness, na.rm=TRUE))
trial_df %>% dplyr::group_by(groupLeth) %>% dplyr::summarize(mean=mean(neo_openness, na.rm=TRUE))

#pid5
trial_df %>% dplyr::group_by(groupLeth) %>% dplyr::summarize(mean=mean(pid5_negativeaffect, na.rm=TRUE))
trial_df %>% dplyr::group_by(groupLeth) %>% dplyr::summarize(mean=mean(pid5_detachment, na.rm=TRUE))
trial_df %>% dplyr::group_by(groupLeth) %>% dplyr::summarize(mean=mean(pid5_antagonism, na.rm=TRUE))
trial_df %>% dplyr::group_by(groupLeth) %>% dplyr::summarize(mean=mean(pid5_disinhibition, na.rm=TRUE))
trial_df %>% dplyr::group_by(groupLeth) %>% dplyr::summarize(mean=mean(pid5_psychoticism, na.rm=TRUE))

anova_groupLeth_pid5depressivity <- aov(Depressivity~group_leth, data=sub_df, na.action=na.exclude)
summary(anova_groupLeth_pid5depressivity)
TukeyHSD(anova_groupLeth_pid5depressivity, na.action=na.exclude)
anova_groupLeth_pid5anhedonia <- aov(Anhedonia~group_leth, data=sub_df, na.action=na.exclude)
summary(anova_groupLeth_pid5anhedonia)
TukeyHSD(anova_groupLeth_pid5anhedonia, na.action=na.exclude)
anova_groupLeth_pid5withdrawal <- aov(Withdrawal~group_leth, data=sub_df, na.action=na.exclude)
summary(anova_groupLeth_pid5withdrawal)
TukeyHSD(anova_groupLeth_pid5withdrawal, na.action=na.exclude)
anova_groupLeth_pid5anx <- aov(Anxiousness~group_leth, data=sub_df, na.action=na.exclude)
summary(anova_groupLeth_pid5anx)
TukeyHSD(anova_groupLeth_pid5anx, na.action=na.exclude)
anova_groupLeth_pid5intavoid <- aov(Intimacy_Avoidance~group_leth, data=sub_df, na.action=na.exclude)
summary(anova_groupLeth_pid5intavoid)
TukeyHSD(anova_groupLeth_pid5intavoid, na.action=na.exclude)
anova_groupLeth_pid5emotlabil <- aov(Emotional_Lability~group_leth, data=sub_df, na.action=na.exclude)
summary(anova_groupLeth_pid5emotlabil)
TukeyHSD(anova_groupLeth_pid5emotlabil, na.action=na.exclude)
anova_groupLeth_pid5persev <- aov(Perseveration~group_leth, data=sub_df, na.action=na.exclude)
summary(anova_groupLeth_pid5persev)
TukeyHSD(anova_groupLeth_pid5persev, na.action=na.exclude)
anova_groupLeth_pid5separinsecur <- aov(Separation_Insecurity~group_leth, data=sub_df, na.action=na.exclude)
summary(anova_groupLeth_pid5separinsecur)
TukeyHSD(anova_groupLeth_pid5separinsecur, na.action=na.exclude)
anova_groupLeth_pid5submis <- aov(Submissiveness~group_leth, data=sub_df, na.action=na.exclude)
summary(anova_groupLeth_pid5submis)
TukeyHSD(anova_groupLeth_pid5submis, na.action=na.exclude)
anova_groupLeth_pid5suspic <- aov(Suspiciousness~group_leth, data=sub_df, na.action=na.exclude)
summary(anova_groupLeth_pid5suspic)
TukeyHSD(anova_groupLeth_pid5suspic, na.action=na.exclude)
anova_groupLeth_pid5restaff <- aov(Restricted_Affectivity~group_leth, data=sub_df, na.action=na.exclude)
summary(anova_groupLeth_pid5restaff)
TukeyHSD(anova_groupLeth_pid5restaff, na.action=na.exclude)
anova_groupLeth_pid5rigid <- aov(Rigid_Perfectionism~group_leth, data=sub_df, na.action=na.exclude)
summary(anova_groupLeth_pid5rigid)
TukeyHSD(anova_groupLeth_pid5rigid, na.action=na.exclude)
anova_groupLeth_pid5grand <- aov(Grandiosity~group_leth, data=sub_df, na.action=na.exclude)
summary(anova_groupLeth_pid5grand)
TukeyHSD(anova_groupLeth_pid5grand, na.action=na.exclude)
anova_groupLeth_pid5risktaking <- aov(Risk_Taking~group_leth, data=sub_df, na.action=na.exclude)
summary(anova_groupLeth_pid5risktaking)
TukeyHSD(anova_groupLeth_pid5risktaking, na.action=na.exclude)
anova_groupLeth_pid5manip <- aov(Manipulativeness~group_leth, data=sub_df, na.action=na.exclude)
summary(anova_groupLeth_pid5manip)
TukeyHSD(anova_groupLeth_pid5manip, na.action=na.exclude)
anova_groupLeth_pid5irrespon <- aov(Irresponsibility~group_leth, data=sub_df, na.action=na.exclude)
summary(anova_groupLeth_pid5irrespon)
TukeyHSD(anova_groupLeth_pid5irrespon, na.action=na.exclude)
anova_groupLeth_pid5impulsivity <- aov(Impulsivity~group_leth, data=sub_df, na.action=na.exclude)
summary(anova_groupLeth_pid5impulsivity)
TukeyHSD(anova_groupLeth_pid5impulsivity, na.action=na.exclude)
anova_groupLeth_pid5hostility <- aov(Hostility~group_leth, data=sub_df, na.action=na.exclude)
summary(anova_groupLeth_pid5hostility)
TukeyHSD(anova_groupLeth_pid5hostility, na.action=na.exclude)
anova_groupLeth_pid5deceitfulness <- aov(Deceitfulness~group_leth, data=sub_df, na.action=na.exclude)
summary(anova_groupLeth_pid5deceitfulness)
TukeyHSD(anova_groupLeth_pid5deceitfulness, na.action=na.exclude)
anova_groupLeth_pid5distract <- aov(Distractibility~group_leth, data=sub_df, na.action=na.exclude)
summary(anova_groupLeth_pid5distract)
TukeyHSD(anova_groupLeth_pid5distract, na.action=na.exclude)
anova_groupLeth_pid5callousness <- aov(Callousness~group_leth, data=sub_df, na.action=na.exclude)
summary(anova_groupLeth_pid5callousness)
TukeyHSD(anova_groupLeth_pid5callousness, na.action=na.exclude)
anova_groupLeth_pid5attentionseeking <- aov(Attention_Seeking~group_leth, data=sub_df, na.action=na.exclude)
summary(anova_groupLeth_pid5attentionseeking)
TukeyHSD(anova_groupLeth_pid5attentionseeking, na.action=na.exclude)
anova_groupLeth_pid5perceptualdysregulation <- aov(Perceptual_Dysregulation~group_leth, data=sub_df, na.action=na.exclude)
summary(anova_groupLeth_pid5perceptualdysregulation)
TukeyHSD(anova_groupLeth_pid5perceptualdysregulation, na.action=na.exclude)
anova_groupLeth_pid5unusualbeliefs <- aov(Unusual_Beliefs_Experiences~group_leth, data=sub_df, na.action=na.exclude)
summary(anova_groupLeth_pid5unusualbeliefs)
TukeyHSD(anova_groupLeth_pid5unusualbeliefs, na.action=na.exclude)
anova_groupLeth_pid5eccentricity <- aov(Eccentricity~group_leth, data=sub_df, na.action=na.exclude)
summary(anova_groupLeth_pid5eccentricity)
TukeyHSD(anova_groupLeth_pid5eccentricity, na.action=na.exclude)
anova_groupLeth_wtar <- aov(wtar~group_leth, data=sub_df, na.action=na.exclude)
summary(anova_groupLeth_wtar)
TukeyHSD(anova_groupLeth_wtar, na.action=na.exclude)
anova_groupLeth_exit <- aov(exit~group_leth, data=sub_df, na.action=na.exclude)
summary(anova_groupLeth_exit)
TukeyHSD(anova_groupLeth_exit, na.action=na.exclude)




###PERSONALITY MEAURES
#Aliona's model
#Lmer (rt ~ rt_lag_sc*reward_lag*scale(TRAIT) + rt_vmax_lag_sc* trial_neg_inv_sc*scale(TRAIT) + rt_lag_sc*trial_neg_inv_sc*scale(TRAIT)+ (1|ID/run), DATASET)

### NEO
#Conscientiousness
#original model
#m10 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(neo_conscientiousness) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(neo_conscientiousness) +
#               rt_lag_sc*condition_trial_neg_inv_sc*scale(neo_conscientiousness) + (1|id/run), trial_df)
#Add random slopes (rt_lag_sc + reward_lag_rec|ID)
#m10 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(neo_conscientiousness) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(neo_conscientiousness) +
#                rt_lag_sc*condition_trial_neg_inv_sc*scale(neo_conscientiousness) + (rt_lag_sc + reward_lag_rec|id) + (1|id/run), trial_df)
#add reward*rt_lag|ID
#m10 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(neo_conscientiousness) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(neo_conscientiousness) +
#                rt_lag_sc*condition_trial_neg_inv_sc*scale(neo_conscientiousness) + (rt_lag_sc*reward_lag_rec|id) + (1|id/run), trial_df)
#Add condition*WSLS 3-way and condition*trait 2-way interactions
#m10 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(neo_conscientiousness) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(neo_conscientiousness) +
#                rt_lag_sc*condition_trial_neg_inv_sc*scale(neo_conscientiousness) + rewFunc*rt_lag_sc*reward_lag_rec +
#                rewFunc*scale(neo_conscientiousness) + (1|id/run), trial_df)

car::Anova(m10, '3')
summary(m10)
avg_trait <- mean(trial_df$neo_conscientiousness, na.rm=TRUE)
sd_trait <- sd(trial_df$neo_conscientiousness, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em10a <- as_tibble(emmeans::emtrends(m10, var = "rt_lag_sc", specs = c("reward_lag_rec", "neo_conscientiousness"), at = list(neo_conscientiousness = c(low_trait, high_trait))))
ggplot(em10a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=neo_conscientiousness)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) +  
    scale_y_reverse()
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em10b <- as_tibble(emmeans::emtrends(m10, var = "rt_vmax_lag_sc", specs = c("neo_conscientiousness", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), neo_conscientiousness = c(low_trait, high_trait)), options = list() ))
ggplot(em10b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=neo_conscientiousness)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))

#Neuroticism
#original model
m11 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(neo_neuroticism) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(neo_neuroticism) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(neo_neuroticism) + (1|id/run), trial_df)
#Add random slope
#m11 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(neo_neuroticism) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(neo_neuroticism) +
#                rt_lag_sc*condition_trial_neg_inv_sc*scale(neo_neuroticism) + (rt_lag_sc + reward_lag_rec|id) + (1|id/run), trial_df)
#add reward*rt_lag|ID
m11 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(neo_neuroticism) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(neo_neuroticism) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(neo_neuroticism) + (rt_lag_sc*reward_lag_rec|id) + (1|id/run), trial_df)
#Add condition*WSLS 3-way and condition*trait 2-way interactions
#m11 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(neo_neuroticism) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(neo_neuroticism) +
#                rt_lag_sc*condition_trial_neg_inv_sc*scale(neo_neuroticism) + rewFunc*rt_lag_sc*reward_lag_rec +
#                rewFunc*scale(neo_neuroticism)  + (1|id/run), trial_df)
car::Anova(m11, '3')
summary(m11)
avg_trait <- mean(trial_df$neo_neuroticism, na.rm=TRUE)
sd_trait <- sd(trial_df$neo_neuroticism, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em11a <- as_tibble(emmeans::emtrends(m11, var = "rt_lag_sc", specs = c("reward_lag_rec", "neo_neuroticism"), at = list(neo_neuroticism = c(low_trait, high_trait))))
ggplot(em11a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=neo_neuroticism)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) +  
    scale_y_reverse()
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em11b <- as_tibble(emmeans::emtrends(m11, var = "rt_vmax_lag_sc", specs = c("neo_neuroticism", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), neo_neuroticism = c(low_trait, high_trait)), options = list() ))
ggplot(em11b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=neo_neuroticism)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))

#Extraversion
#original model
#m12 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(neo_extraversion) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(neo_extraversion) +
#                rt_lag_sc*condition_trial_neg_inv_sc*scale(neo_extraversion) + (1|id/run), trial_df)
#Add random slope
#m12 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(neo_extraversion) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(neo_extraversion) +
#                rt_lag_sc*condition_trial_neg_inv_sc*scale(neo_extraversion) + (rt_lag_sc + reward_lag_rec|id) + (1|id/run), trial_df)
#add reward*rt_lag|ID
m12 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(neo_extraversion) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(neo_extraversion) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(neo_extraversion) + (rt_lag_sc*reward_lag_rec|id) + (1|id/run), trial_df)
#Add condition*WSLS 3-way and condition*trait 2-way interactions
#m12 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(neo_extraversion) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(neo_extraversion) +
#                rt_lag_sc*condition_trial_neg_inv_sc*scale(neo_extraversion) + rewFunc*rt_lag_sc*reward_lag_rec +
#                rewFunc*scale(neo_extraversion) + (1|id/run), trial_df)
car::Anova(m12, '3')
summary(m12)
avg_trait <- mean(trial_df$neo_extraversion, na.rm=TRUE)
sd_trait <- sd(trial_df$neo_extraversion, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em12a <- as_tibble(emmeans::emtrends(m12, var = "rt_lag_sc", specs = c("reward_lag_rec", "neo_extraversion"), at = list(neo_extraversion = c(low_trait, high_trait))))
ggplot(em12a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=neo_extraversion)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) +  
    scale_y_reverse()
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em12b <- as_tibble(emmeans::emtrends(m12, var = "rt_vmax_lag_sc", specs = c("neo_extraversion", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), neo_extraversion = c(low_trait, high_trait)), options = list() ))
ggplot(em12b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=neo_extraversion)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))

#Agreeableness
#original model
#m13 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(neo_aggreableness) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(neo_aggreableness) +
#                rt_lag_sc*condition_trial_neg_inv_sc*scale(neo_aggreableness) + (1|id/run), trial_df)
#Add random slopes
#m13 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(neo_aggreableness) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(neo_aggreableness) +
#                rt_lag_sc*condition_trial_neg_inv_sc*scale(neo_aggreableness) + (rt_lag_sc + reward_lag_rec|id) + (1|id/run), trial_df)
#add reward*rt_lag|ID
m13 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(neo_aggreableness) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(neo_aggreableness) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(neo_aggreableness) + (rt_lag_sc*reward_lag_rec|id) + (1|id/run), trial_df)
#Add condition*WSLS 3-way and condition*trait 2-way interactions
#m13 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(neo_aggreableness) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(neo_aggreableness) +
#                rt_lag_sc*condition_trial_neg_inv_sc*scale(neo_aggreableness) + rewFunc*rt_lag_sc*reward_lag_rec +
#                rewFunc*scale(neo_aggreableness) + (1|id/run), trial_df)
car::Anova(m13, '3')
summary(m13)
avg_trait <- mean(trial_df$neo_aggreableness, na.rm=TRUE)
sd_trait <- sd(trial_df$neo_aggreableness, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em13a <- as_tibble(emmeans::emtrends(m13, var = "rt_lag_sc", specs = c("reward_lag_rec", "neo_aggreableness"), at = list(neo_aggreableness = c(low_trait, high_trait))))
ggplot(em13a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=neo_aggreableness)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) +  
    scale_y_reverse()
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em13b <- as_tibble(emmeans::emtrends(m13, var = "rt_vmax_lag_sc", specs = c("neo_aggreableness", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), neo_aggreableness = c(low_trait, high_trait)), options = list() ))
ggplot(em13b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=neo_aggreableness)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))

#Openness
#original model
#m14 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(neo_openness) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(neo_openness) +
#                rt_lag_sc*condition_trial_neg_inv_sc*scale(neo_openness) + (1|id/run), trial_df)
#with random slopes
#m14 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(neo_openness) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(neo_openness) +
#                rt_lag_sc*condition_trial_neg_inv_sc*scale(neo_openness) + (rt_lag_sc + reward_lag_rec|id) + (1|id/run), trial_df)
#add reward*rt_lag|ID
m14 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(neo_openness) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(neo_openness) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(neo_openness) + (rt_lag_sc*reward_lag_rec|id) + (1|id/run), trial_df)
#Add condition*WSLS 3-way and condition*trait 2-way interactions
#m14 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(neo_openness) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(neo_openness) +
#                rt_lag_sc*condition_trial_neg_inv_sc*scale(neo_openness) + rewFunc*rt_lag_sc*reward_lag_rec +
#                rewFunc*scale(neo_openness) + (1|id/run), trial_df)
car::Anova(m14, '3')
summary(m14)
avg_trait <- mean(trial_df$neo_openness, na.rm=TRUE)
sd_trait <- sd(trial_df$neo_openness, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em14a <- as_tibble(emmeans::emtrends(m14, var = "rt_lag_sc", specs = c("reward_lag_rec", "neo_openness"), at = list(neo_openness = c(low_trait, high_trait))))
ggplot(em14a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=neo_openness)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) +  
    scale_y_reverse()
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em14b <- as_tibble(emmeans::emtrends(m14, var = "rt_vmax_lag_sc", specs = c("neo_openness", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), neo_openness = c(low_trait, high_trait)), options = list() ))
ggplot(em14b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=neo_openness)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))

###UPPSP Personality Measures
#total
m15_tot <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(uppsp_total) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(uppsp_total) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(uppsp_total) + (1|id/run), trial_df)
car::Anova(m15_tot, '3')
summary(m15_tot)
avg_trait <- mean(trial_df$uppsp_total, na.rm=TRUE)
sd_trait <- sd(trial_df$uppsp_total, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em15_tot <- as_tibble(emmeans::emtrends(m15_tot, var = "rt_lag_sc", specs = c("reward_lag_rec", "uppsp_total"), at = list(uppsp_total = c(low_trait, high_trait))))
ggplot(em15_tot, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=uppsp_total)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) +  
    scale_y_reverse()
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em15_totb <- as_tibble(emmeans::emtrends(m15_tot, var = "rt_vmax_lag_sc", specs = c("uppsp_total", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), uppsp_total = c(low_trait, high_trait)), options = list() ))
ggplot(em15_totb, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=uppsp_total)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))

### October 2023 - Extract coefficients - to use in brain to behavior analysis
#Maybe want to create new model to independently extract how much people tend to shift after losses, rewards and how much people converge on best RT early vs. late
trial_df_loss <- trial_df %>% filter(reward_lag_rec=="-0.5")
trial_df_win <- trial_df %>% filter(reward_lag_rec=="0.5")
trial_df_late <- trial_df %>% filter(trial_df$condition_trial>20)

trial_df_corr <- trial_df %>% group_by(id) %>%
    mutate(pearson_rt_corr = cor(rt_csv,rt_lag_sc, method = "pearson",use="complete.obs"),
           pearson_rtvmax_corr = cor(rt_csv,rt_vmax_lag_sc, method = "pearson", use="complete.obs"))
trial_df_corr_late <- trial_df_late %>% group_by(id) %>%
    mutate(pearson_rt_corr = cor(rt_csv,rt_lag_sc, method = "pearson",use="complete.obs"),
           pearson_rtvmax_corr = cor(rt_csv,rt_vmax_lag_sc, method = "pearson", use="complete.obs"))
trial_df_corr_loss <- trial_df_loss %>% group_by(id) %>%
    mutate(pearson_rt_corr = cor(rt_csv,rt_lag_sc, method = "pearson",use="complete.obs"),
           pearson_rtvmax_corr = cor(rt_csv,rt_vmax_lag_sc, method = "pearson",use="complete.obs"))
trial_df_corr_win <- trial_df_win %>% group_by(id) %>%
    mutate(pearson_rt_corr = cor(rt_csv,rt_lag_sc, method = "pearson",use="complete.obs"),
           pearson_rtvmax_corr = cor(rt_csv,rt_vmax_lag_sc, method = "pearson",use="complete.obs"))

#Put these new regressors into demos dataframe (Oct 19, 2023)
demos_corr_whole <- trial_df_corr %>% mutate(pearson_rtvmax_corr_whole = pearson_rtvmax_corr,
                                             pearson_rt_corr_whole = pearson_rt_corr) %>% 
    select(id, pearson_rt_corr_whole, pearson_rtvmax_corr_whole) %>% unique()
demos_corr_late <- trial_df_corr_late  %>% mutate(pearson_rtvmax_corr_late = pearson_rtvmax_corr,
                                                  pearson_rt_corr_late = pearson_rt_corr) %>% 
    select(id, pearson_rt_corr_late, pearson_rtvmax_corr_late) %>% unique()
demos_corr_loss <- trial_df_corr_loss %>% mutate(pearson_rtvmax_corr_loss = pearson_rtvmax_corr,
                                                 pearson_rt_corr_loss = pearson_rt_corr) %>% 
    select(id, pearson_rt_corr_loss, pearson_rtvmax_corr_loss) %>% unique() 
demos_corr_win <- trial_df_corr_win %>% mutate(pearson_rtvmax_corr_win = pearson_rtvmax_corr,
                                               pearson_rt_corr_win = pearson_rt_corr) %>% 
    select(id, pearson_rt_corr_win, pearson_rtvmax_corr_win) %>% unique() 

demos_corr_compiled <- demos_corr_whole %>% left_join(demos_corr_late, by="id")
demos_corr_compiled <- demos_corr_compiled %>% left_join(demos_corr_loss, by="id")
demos_corr_compiled <- demos_corr_compiled %>% left_join(demos_corr_win, by="id")
demos_corr_compiled$id <- as.character(demos_corr_compiled$id)
demos <- demos %>% left_join(demos_corr_compiled, by="id")


demos %>% group_by(groupLeth) %>% summarise(mean_rt_corr_whole=mean(pearson_rt_corr_whole),
                                            mean_rt_corr_win=mean(pearson_rt_corr_win),
                                            mean_rt_corr_loss=mean(pearson_rt_corr_loss),
                                            mean_rtvmax_corr_whole=mean(pearson_rtvmax_corr_whole),
                                            mean_rtvmax_corr_late=mean(pearson_rtvmax_corr_late))

##RT swings
#Impulsivity corr - simple_groups, Group, groupLeth
trial_df_corr %>% group_by(simple_groups) %>% select(id, groupLeth, pearson_rt_corr, uppsp_total, neo_neuroticism) %>% 
    unique() %>% mutate(pearson_rt_corr_neg = -pearson_rt_corr) %>%
    summarise(impulsivity_rtswing_corr = cor(uppsp_total,pearson_rt_corr_neg,use="complete.obs"))
trial_df_corr_late %>% group_by(simple_groups) %>% select(id, groupLeth, pearson_rt_corr, uppsp_total, neo_neuroticism) %>% 
    unique() %>% mutate(pearson_rt_corr_neg = -pearson_rt_corr) %>%
    summarise(impulsivity_rtswing_corr = cor(uppsp_total,pearson_rt_corr_neg,use="complete.obs"))
trial_df_corr_win %>% group_by(simple_groups) %>% select(id, groupLeth, pearson_rt_corr, uppsp_total, neo_neuroticism) %>% 
    unique() %>% mutate(pearson_rt_corr_neg = -pearson_rt_corr) %>%
    summarise(impulsivity_rtswing_corr = cor(uppsp_total,pearson_rt_corr_neg,use="complete.obs"))
trial_df_corr_loss %>% group_by(simple_groups) %>% select(id, groupLeth, pearson_rt_corr, uppsp_total, neo_neuroticism) %>% 
    unique() %>% mutate(pearson_rt_corr_neg = -pearson_rt_corr) %>%
    summarise(impulsivity_rtswing_corr = cor(uppsp_total,pearson_rt_corr_neg,use="complete.obs"))

#Neuroticism corr
trial_df_corr %>% group_by(simple_groups) %>% select(id, groupLeth, pearson_rt_corr, uppsp_total, neo_neuroticism) %>% 
    unique() %>% mutate(pearson_rt_corr_neg = -pearson_rt_corr) %>%
    summarise(neuroticism_rtswing_corr = cor(neo_neuroticism,pearson_rt_corr_neg,use="complete.obs"))
trial_df_corr_late %>% group_by(simple_groups) %>% select(id, groupLeth, pearson_rt_corr, uppsp_total, neo_neuroticism) %>% 
    unique() %>% mutate(pearson_rt_corr_neg = -pearson_rt_corr) %>%
    summarise(neuroticism_rtswing_corr = cor(neo_neuroticism,pearson_rt_corr_neg,use="complete.obs"))
trial_df_corr_win %>% group_by(simple_groups) %>% select(id, groupLeth, pearson_rt_corr, uppsp_total, neo_neuroticism) %>% 
    unique() %>% mutate(pearson_rt_corr_neg = -pearson_rt_corr) %>%
    summarise(neuroticism_rtswing_corr = cor(neo_neuroticism,pearson_rt_corr_neg,use="complete.obs"))
trial_df_corr_loss %>% group_by(simple_groups) %>% select(id, groupLeth, pearson_rt_corr, uppsp_total, neo_neuroticism) %>% 
    unique() %>% mutate(pearson_rt_corr_neg = -pearson_rt_corr) %>%
    summarise(neuroticism_rtswing_corr = cor(neo_neuroticism,pearson_rt_corr_neg,use="complete.obs"))

##Exploitation
#Impulsivity corr - simple_groups, Group, groupLeth
trial_df_corr %>% group_by(simple_groups) %>% select(id, groupLeth, pearson_rtvmax_corr, uppsp_total, neo_neuroticism) %>% 
    unique() %>% 
    summarise(impulsivity_rtvmax_corr = cor(uppsp_total, pearson_rtvmax_corr,use="complete.obs"))
trial_df_corr_late %>% group_by(simple_groups) %>% select(id, groupLeth, pearson_rtvmax_corr, uppsp_total, neo_neuroticism) %>% 
    unique() %>% 
    summarise(impulsivity_rtvmax_corr = cor(uppsp_total, pearson_rtvmax_corr,use="complete.obs"))
trial_df_corr_win %>% group_by(simple_groups) %>% select(id, groupLeth, pearson_rtvmax_corr, uppsp_total, neo_neuroticism) %>% 
    unique() %>% 
    summarise(impulsivity_rtvmax_corr = cor(uppsp_total,pearson_rtvmax_corr,use="complete.obs"))
trial_df_corr_loss %>% group_by(simple_groups) %>% select(id, groupLeth, pearson_rtvmax_corr, uppsp_total, neo_neuroticism) %>% 
    unique() %>% 
    summarise(impulsivity_rtvmax_corr = cor(uppsp_total,pearson_rtvmax_corr,use="complete.obs"))

#Neuroticism corr
trial_df_corr %>% group_by(simple_groups) %>% select(id, groupLeth, pearson_rtvmax_corr, uppsp_total, neo_neuroticism) %>% 
    unique() %>% 
    summarise(neuroticism_rtvmax_corr = cor(neo_neuroticism,pearson_rtvmax_corr,use="complete.obs"))
trial_df_corr_late %>% group_by(simple_groups) %>% select(id, groupLeth, pearson_rtvmax_corr, uppsp_total, neo_neuroticism) %>% 
    unique() %>% 
    summarise(neuroticism_rtvmax_corr = cor(neo_neuroticism,pearson_rtvmax_corr,use="complete.obs"))
trial_df_corr_win %>% group_by(simple_groups) %>% select(id, groupLeth, pearson_rtvmax_corr, uppsp_total, neo_neuroticism) %>% 
    unique() %>% 
    summarise(neuroticism_rtvmax_corr = cor(neo_neuroticism,pearson_rtvmax_corr,use="complete.obs"))
trial_df_corr_loss %>% group_by(simple_groups) %>% select(id, groupLeth, pearson_rtvmax_corr, uppsp_total, neo_neuroticism) %>% 
    unique() %>% 
    summarise(neuroticism_rtvmax_corr = cor(neo_neuroticism,pearson_rtvmax_corr,use="complete.obs"))


#trial_df_corr_loss %>% summarise(impulsivity_rtswing_corr = cor(uppsp_total,pearson_rt_corr,use="complete.obs"))
#rt_loss.aov <- aov(pearson_rt_corr ~ groupLeth, data = trial_df_corr_loss)
#summary(rt_loss.aov)
#Anova(trial_df_corr, '3')
#summary(m15)

#summarise(pearson = cor(rt_csv,rt_lag_sc, method = "pearson"),
#          spearman = cor(rt_csv,rt_lag_sc, method = "spearman")) %>%
#    filter(reward_lag_rec==0.5)

#Doesn't work:
#m100_loss <- lmer(rt_csv ~ rt_lag_sc + rt_vmax_lag_sc*condition_trial_neg_inv_sc +
#                      rt_lag_sc*condition_trial_neg_inv_sc + (1|id/run), trial_df_loss)
#m100_win <- lmer(rt_csv ~ rt_lag_sc + rt_vmax_lag_sc*condition_trial_neg_inv_sc +
#                     rt_lag_sc*condition_trial_neg_inv_sc + (1|id/run), trial_df_win)
#loss_coefs <- coef(m100_loss)
#win_coefs <- coef(m100_win)

#Negative Urgency
m15 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(uppsp_negative_urgency) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(uppsp_negative_urgency) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(uppsp_negative_urgency) + (1|id/run), trial_df)
Anova(m15, '3')
summary(m15)

#Positive Urgency
m16 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(uppsp_positive_urgency) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(uppsp_positive_urgency) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(uppsp_positive_urgency) + (1|id/run), trial_df)
Anova(m16, '3')
summary(m16)

#Lack of Premeditation
m17 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(uppsp_lack_of_premeditation) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(uppsp_lack_of_premeditation) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(uppsp_lack_of_premeditation) + (1|id/run), trial_df)
Anova(m17, '3')
summary(m17)
avg_trait <- mean(trial_df$uppsp_lack_of_premeditation, na.rm=TRUE)
sd_trait <- sd(trial_df$uppsp_lack_of_premeditation, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em17a <- as_tibble(emmeans::emtrends(m17, var = "rt_lag_sc", specs = c("reward_lag_rec", "uppsp_lack_of_premeditation"), at = list(uppsp_lack_of_premeditation = c(low_trait, high_trait))))
ggplot(em17a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=uppsp_lack_of_premeditation)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) +  
    scale_y_reverse()
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em17b <- as_tibble(emmeans::emtrends(m17, var = "rt_vmax_lag_sc", specs = c("uppsp_lack_of_premeditation", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), uppsp_lack_of_premeditation = c(low_trait, high_trait)), options = list() ))
ggplot(em17b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=uppsp_lack_of_premeditation)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))

#Lack of Perseveration
m18 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(uppsp_lack_of_perseveration) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(uppsp_lack_of_perseveration) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(uppsp_lack_of_perseveration) + (1|id/run), trial_df)
Anova(m18, '3')
summary(m18)
avg_trait <- mean(trial_df$uppsp_lack_of_perseveration, na.rm=TRUE)
sd_trait <- sd(trial_df$uppsp_lack_of_perseveration, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em18a <- as_tibble(emmeans::emtrends(m18, var = "rt_lag_sc", specs = c("reward_lag_rec", "uppsp_lack_of_perseveration"), at = list(uppsp_lack_of_perseveration = c(low_trait, high_trait))))
ggplot(em18a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=uppsp_lack_of_perseveration)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) +  
    scale_y_reverse()
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em18b <- as_tibble(emmeans::emtrends(m18, var = "rt_vmax_lag_sc", specs = c("uppsp_lack_of_perseveration", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), uppsp_lack_of_perseveration = c(low_trait, high_trait)), options = list() ))
ggplot(em18b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=uppsp_lack_of_perseveration)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))

###PID5
##Internalizing
#Depressivity, pid5_depressivity
m20 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_depressivity) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_depressivity) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_depressivity) + (1|id/run), trial_df)
Anova(m20, '3')
summary(m20)
avg_trait <- mean(trial_df$pid5_depressivity, na.rm=TRUE)
sd_trait <- sd(trial_df$pid5_depressivity, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em20a <- as_tibble(emmeans::emtrends(m20, var = "rt_lag_sc", specs = c("reward_lag_rec", "pid5_depressivity"), at = list(pid5_depressivity = c(low_trait, high_trait))))
ggplot(em20a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_depressivity)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) +  
    scale_y_reverse()
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em20b <- as_tibble(emmeans::emtrends(m20, var = "rt_vmax_lag_sc", specs = c("pid5_depressivity", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), pid5_depressivity = c(low_trait, high_trait)), options = list() ))
ggplot(em20b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_depressivity)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))

#Anhedonia, pid5_anhedonia
m21 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_anhedonia) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_anhedonia) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_anhedonia) + (1|id/run), trial_df)
Anova(m21, '3')
summary(m21)
avg_trait <- mean(trial_df$pid5_anhedonia, na.rm=TRUE)
sd_trait <- sd(trial_df$pid5_anhedonia, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em21a <- as_tibble(emmeans::emtrends(m21, var = "rt_lag_sc", specs = c("reward_lag_rec", "pid5_anhedonia"), at = list(pid5_anhedonia = c(low_trait, high_trait))))
ggplot(em21a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_anhedonia)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) +  
    scale_y_reverse()
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em21b <- as_tibble(emmeans::emtrends(m21, var = "rt_vmax_lag_sc", specs = c("pid5_anhedonia", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), pid5_anhedonia = c(low_trait, high_trait)), options = list() ))
ggplot(em21b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_anhedonia)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))

#Withdrawal, pid5_withdrawal
m22 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_withdrawal) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_withdrawal) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_withdrawal) + (1|id/run), trial_df)
Anova(m22, '3')
summary(m22)
avg_trait <- mean(trial_df$pid5_withdrawal, na.rm=TRUE)
sd_trait <- sd(trial_df$pid5_withdrawal, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em22a <- as_tibble(emmeans::emtrends(m22, var = "rt_lag_sc", specs = c("reward_lag_rec", "pid5_withdrawal"), at = list(pid5_withdrawal = c(low_trait, high_trait))))
ggplot(em22a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_withdrawal)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) +  
    scale_y_reverse()
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em22b <- as_tibble(emmeans::emtrends(m22, var = "rt_vmax_lag_sc", specs = c("pid5_withdrawal", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), pid5_withdrawal = c(low_trait, high_trait)), options = list() ))
ggplot(em22b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_withdrawal)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))

#Anxiousness, pid5_anxiousness
m23 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_anxiousness) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_anxiousness) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_anxiousness) + (1|id/run), trial_df)
Anova(m23, '3')
summary(m23)
avg_trait <- mean(trial_df$pid5_anxiousness, na.rm=TRUE)
sd_trait <- sd(trial_df$pid5_anxiousness, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em23a <- as_tibble(emmeans::emtrends(m23, var = "rt_lag_sc", specs = c("reward_lag_rec", "pid5_anxiousness"), at = list(pid5_anxiousness = c(low_trait, high_trait))))
ggplot(em23a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_anxiousness)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) +  
    scale_y_reverse()
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em23b <- as_tibble(emmeans::emtrends(m23, var = "rt_vmax_lag_sc", specs = c("pid5_anxiousness", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), pid5_anxiousness = c(low_trait, high_trait)), options = list() ))
ggplot(em23b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_anxiousness)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))

#Intimacy Avoidance, pid5_intimacy_avoidance - n.s.
m24 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_intimacy_avoidance) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_intimacy_avoidance) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_intimacy_avoidance) + (1|id/run), trial_df)
Anova(m24, '3')

#Emotional Lability, pid5_emotional_lability
m25 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_emotional_lability) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_emotional_lability) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_emotional_lability) + (1|id/run), trial_df)
Anova(m25, '3')
summary(m25)
avg_trait <- mean(trial_df$pid5_emotional_lability, na.rm=TRUE)
sd_trait <- sd(trial_df$pid5_emotional_lability, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em25a <- as_tibble(emmeans::emtrends(m25, var = "rt_lag_sc", specs = c("reward_lag_rec", "pid5_emotional_lability"), at = list(pid5_emotional_lability = c(low_trait, high_trait))))
ggplot(em25a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_emotional_lability)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) +  
    scale_y_reverse()
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em25b <- as_tibble(emmeans::emtrends(m25, var = "rt_vmax_lag_sc", specs = c("pid5_emotional_lability", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), pid5_emotional_lability = c(low_trait, high_trait)), options = list() ))
ggplot(em25b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_emotional_lability)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))

#Perseveration, pid5_perseveration
m26 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_perseveration) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_perseveration) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_perseveration) + (1|id/run), trial_df)
Anova(m26, '3')
summary(m26)
avg_trait <- mean(trial_df$pid5_perseveration, na.rm=TRUE)
sd_trait <- sd(trial_df$pid5_perseveration, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em26a <- as_tibble(emmeans::emtrends(m26, var = "rt_lag_sc", specs = c("reward_lag_rec", "pid5_perseveration"), at = list(pid5_perseveration = c(low_trait, high_trait))))
ggplot(em26a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_perseveration)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) +  
    scale_y_reverse()
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em26b <- as_tibble(emmeans::emtrends(m26, var = "rt_vmax_lag_sc", specs = c("pid5_perseveration", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), pid5_perseveration = c(low_trait, high_trait)), options = list() ))
ggplot(em26b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_perseveration)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))

#Separation insecurity, pid5_separation_insecurity - n.s.
m27 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_separation_insecurity) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_separation_insecurity) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_separation_insecurity) + (1|id/run), trial_df)
Anova(m27, '3')

#Submissiveness, pid5_submissiveness - n.s.
m28 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_submissiveness) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_submissiveness) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_submissiveness) + (1|id/run), trial_df)
Anova(m28, '3')

#Suspiciousness, pid5_suspiciousness
m29 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_suspiciousness) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_suspiciousness) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_suspiciousness) + (1|id/run), trial_df)
Anova(m29, '3')
summary(m29)
avg_trait <- mean(trial_df$pid5_suspiciousness, na.rm=TRUE)
sd_trait <- sd(trial_df$pid5_suspiciousness, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em29a <- as_tibble(emmeans::emtrends(m29, var = "rt_lag_sc", specs = c("reward_lag_rec", "pid5_suspiciousness"), at = list(pid5_suspiciousness = c(low_trait, high_trait))))
ggplot(em29a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_suspiciousness)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) +  
    scale_y_reverse()
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em29b <- as_tibble(emmeans::emtrends(m29, var = "rt_vmax_lag_sc", specs = c("pid5_suspiciousness", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), pid5_suspiciousness = c(low_trait, high_trait)), options = list() ))
ggplot(em29b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_suspiciousness)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))

#Restricted affectivity, pid5_restricted_affectivity - n.s.
m30 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_restricted_affectivity) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_restricted_affectivity) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_restricted_affectivity) + (1|id/run), trial_df)
Anova(m30, '3')

#Rigid perfectionism, pid5_rigid_perfectionism
m31 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_rigid_perfectionism) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_rigid_perfectionism) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_rigid_perfectionism) + (1|id/run), trial_df)
Anova(m31, '3')
summary(m31)
avg_trait <- mean(trial_df$pid5_rigid_perfectionism, na.rm=TRUE)
sd_trait <- sd(trial_df$pid5_rigid_perfectionism, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em31a <- as_tibble(emmeans::emtrends(m31, var = "rt_lag_sc", specs = c("reward_lag_rec", "pid5_rigid_perfectionism"), at = list(pid5_rigid_perfectionism = c(low_trait, high_trait))))
ggplot(em31a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_rigid_perfectionism)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) +  
    scale_y_reverse()
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em31b <- as_tibble(emmeans::emtrends(m31, var = "rt_vmax_lag_sc", specs = c("pid5_rigid_perfectionism", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), pid5_rigid_perfectionism = c(low_trait, high_trait)), options = list() ))
ggplot(em31b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_rigid_perfectionism)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))

##Externalizing 
#Grandiosity, pid5_grandiosity
m32 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_grandiosity) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_grandiosity) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_grandiosity) + (1|id/run), trial_df)
Anova(m32, '3')
summary(m32)
avg_trait <- mean(trial_df$pid5_grandiosity, na.rm=TRUE)
sd_trait <- sd(trial_df$pid5_grandiosity, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em32a <- as_tibble(emmeans::emtrends(m32, var = "rt_lag_sc", specs = c("reward_lag_rec", "pid5_grandiosity"), at = list(pid5_grandiosity = c(low_trait, high_trait))))
ggplot(em32a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_grandiosity)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) +  
    scale_y_reverse()
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em32b <- as_tibble(emmeans::emtrends(m32, var = "rt_vmax_lag_sc", specs = c("pid5_grandiosity", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), pid5_grandiosity = c(low_trait, high_trait)), options = list() ))
ggplot(em32b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_grandiosity)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))

#Risk taking, pid5_risk_taking
m33 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_risk_taking) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_risk_taking) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_risk_taking) + (1|id/run), trial_df)
Anova(m33, '3')
summary(m33)
avg_trait <- mean(trial_df$pid5_risk_taking, na.rm=TRUE)
sd_trait <- sd(trial_df$pid5_risk_taking, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em33a <- as_tibble(emmeans::emtrends(m33, var = "rt_lag_sc", specs = c("reward_lag_rec", "pid5_risk_taking"), at = list(pid5_risk_taking = c(low_trait, high_trait))))
ggplot(em33a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_risk_taking)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) +  
    scale_y_reverse()
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em33b <- as_tibble(emmeans::emtrends(m33, var = "rt_vmax_lag_sc", specs = c("pid5_risk_taking", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), pid5_risk_taking = c(low_trait, high_trait)), options = list() ))
ggplot(em33b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_risk_taking)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))

#Manipululativeness, pid5_manipulativeness
m34 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_manipulativeness) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_manipulativeness) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_manipulativeness) + (1|id/run), trial_df)
Anova(m34, '3')
summary(m34)
avg_trait <- mean(trial_df$pid5_manipulativeness, na.rm=TRUE)
sd_trait <- sd(trial_df$pid5_manipulativeness, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em34a <- as_tibble(emmeans::emtrends(m34, var = "rt_lag_sc", specs = c("reward_lag_rec", "pid5_manipulativeness"), at = list(pid5_manipulativeness = c(low_trait, high_trait))))
ggplot(em34a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_manipulativeness)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) +  
    scale_y_reverse()
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em34b <- as_tibble(emmeans::emtrends(m34, var = "rt_vmax_lag_sc", specs = c("pid5_manipulativeness", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), pid5_manipulativeness = c(low_trait, high_trait)), options = list() ))
ggplot(em34b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_manipulativeness)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))

#Irresponsibility, pid5_irresponsibility - exploitation only
m35 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_irresponsibility) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_irresponsibility) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_irresponsibility) + (1|id/run), trial_df)
Anova(m35, '3')
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em35b <- as_tibble(emmeans::emtrends(m35, var = "rt_vmax_lag_sc", specs = c("pid5_irresponsibility", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), pid5_irresponsibility = c(low_trait, high_trait)), options = list() ))
ggplot(em35b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_irresponsibility)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))

#Impulsivity, pid5_impulsivity -
m36 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_impulsivity) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_impulsivity) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_impulsivity) + (1|id/run), trial_df)
Anova(m36, '3')

#Hostility, pid5_hostility
m37 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_hostility) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_hostility) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_hostility) + (1|id/run), trial_df)
Anova(m37, '3')
summary(m37)
avg_trait <- mean(trial_df$pid5_hostility, na.rm=TRUE)
sd_trait <- sd(trial_df$pid5_hostility, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em37a <- as_tibble(emmeans::emtrends(m37, var = "rt_lag_sc", specs = c("reward_lag_rec", "pid5_hostility"), at = list(pid5_hostility = c(low_trait, high_trait))))
ggplot(em37a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_hostility)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) +  
    scale_y_reverse()
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em37b <- as_tibble(emmeans::emtrends(m37, var = "rt_vmax_lag_sc", specs = c("pid5_hostility", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), pid5_hostility = c(low_trait, high_trait)), options = list() ))
ggplot(em37b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_hostility)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))

#Deceitfulness, pid5_deceitfulness - exploitation effects only
m38 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_deceitfulness) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_deceitfulness) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_deceitfulness) + (1|id/run), trial_df)
Anova(m38, '3')
#Plot exploitation
avg_trait <- mean(trial_df$pid5_deceitfulness, na.rm=TRUE)
sd_trait <- sd(trial_df$pid5_deceitfulness, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em38b <- as_tibble(emmeans::emtrends(m38, var = "rt_vmax_lag_sc", specs = c("pid5_deceitfulness", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), pid5_deceitfulness = c(low_trait, high_trait)), options = list() ))
ggplot(em38b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_deceitfulness)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))

#Distractability, pid5_distractibility - n.s.
m39 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_distractibility) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_distractibility) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_distractibility) + (1|id/run), trial_df)
Anova(m39, '3')

#Callousness, pid5_callousness
m40 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_callousness) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_callousness) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_callousness) + (1|id/run), trial_df)
Anova(m40, '3')
summary(m40)
avg_trait <- mean(trial_df$pid5_callousness, na.rm=TRUE)
sd_trait <- sd(trial_df$pid5_callousness, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em40a <- as_tibble(emmeans::emtrends(m40, var = "rt_lag_sc", specs = c("reward_lag_rec", "pid5_callousness"), at = list(pid5_callousness = c(low_trait, high_trait))))
ggplot(em40a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_callousness)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) +  
    scale_y_reverse()
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em40b <- as_tibble(emmeans::emtrends(m40, var = "rt_vmax_lag_sc", specs = c("pid5_callousness", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), pid5_callousness = c(low_trait, high_trait)), options = list() ))
ggplot(em40b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_callousness)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))

#Attention Seeking, pid5_attention_seeking
m41 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_attention_seeking) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_attention_seeking) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_attention_seeking) + (1|id/run), trial_df)
Anova(m41, '3')
#Plot exploitation
avg_trait <- mean(trial_df$pid5_attention_seeking, na.rm=TRUE)
sd_trait <- sd(trial_df$pid5_attention_seeking, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em41b <- as_tibble(emmeans::emtrends(m41, var = "rt_vmax_lag_sc", specs = c("pid5_attention_seeking", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), pid5_attention_seeking = c(low_trait, high_trait)), options = list() ))
ggplot(em41b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_attention_seeking)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))


##Psychosis 
#Perceptual dysregulation, pid5_perceptual_dysregulation
m42 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_perceptual_dysregulation) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_perceptual_dysregulation) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_perceptual_dysregulation) + (1|id/run), trial_df)
Anova(m42, '3')
summary(m42)
avg_trait <- mean(trial_df$pid5_perceptual_dysregulation, na.rm=TRUE)
sd_trait <- sd(trial_df$pid5_perceptual_dysregulation, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em42a <- as_tibble(emmeans::emtrends(m42, var = "rt_lag_sc", specs = c("reward_lag_rec", "pid5_perceptual_dysregulation"), at = list(pid5_perceptual_dysregulation = c(low_trait, high_trait))))
ggplot(em42a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_perceptual_dysregulation)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) +  
    scale_y_reverse()
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em42b <- as_tibble(emmeans::emtrends(m42, var = "rt_vmax_lag_sc", specs = c("pid5_perceptual_dysregulation", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), pid5_perceptual_dysregulation = c(low_trait, high_trait)), options = list() ))
ggplot(em42b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_perceptual_dysregulation)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))

#Unusual beliefs, pid5_unusual_beliefs_experiences
m43 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_unusual_beliefs_experiences) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_unusual_beliefs_experiences) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_unusual_beliefs_experiences) + (1|id/run), trial_df)
Anova(m43, '3')
summary(m43)
avg_trait <- mean(trial_df$pid5_unusual_beliefs_experiences, na.rm=TRUE)
sd_trait <- sd(trial_df$pid5_unusual_beliefs_experiences, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em43a <- as_tibble(emmeans::emtrends(m43, var = "rt_lag_sc", specs = c("reward_lag_rec", "pid5_unusual_beliefs_experiences"), at = list(pid5_unusual_beliefs_experiences = c(low_trait, high_trait))))
ggplot(em43a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_unusual_beliefs_experiences)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) +  
    scale_y_reverse()

#Eccentricity, pid5_eccentricity
m44 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_eccentricity) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_eccentricity) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_eccentricity) + (1|id/run), trial_df)
Anova(m44, '3')

### PID5 subcomponent effects
#Negative affect
#original model
#m45 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_negativeaffect) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_negativeaffect) +
#                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_negativeaffect) + (1|id/run), trial_df)
#add random slopes
#m45 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_negativeaffect) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_negativeaffect) +
#                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_negativeaffect) + (rt_lag_sc + reward_lag_rec|id) + (1|id/run), trial_df)
#add reward*rt_lag|ID
m45 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_negativeaffect) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_negativeaffect) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_negativeaffect) + (rt_lag_sc*reward_lag_rec|id) + (1|id/run), trial_df)
#Add condition*WSLS 3-way and condition*trait 2-way interactions
#m45 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_negativeaffect) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_negativeaffect) +
#                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_negativeaffect) + rewFunc*rt_lag_sc*reward_lag_rec +
#                rewFunc*scale(pid5_negativeaffect) + (1|id/run), trial_df)
car::Anova(m45, '3')
summary(m45)
avg_trait <- mean(trial_df$pid5_negativeaffect, na.rm=TRUE)
sd_trait <- sd(trial_df$pid5_negativeaffect, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em45a <- as_tibble(emmeans::emtrends(m45, var = "rt_lag_sc", specs = c("reward_lag_rec", "pid5_negativeaffect"), at = list(pid5_negativeaffect = c(low_trait, high_trait))))
ggplot(em45a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_negativeaffect)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) +  
    scale_y_reverse()
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em45b <- as_tibble(emmeans::emtrends(m45, var = "rt_vmax_lag_sc", specs = c("pid5_negativeaffect", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), pid5_negativeaffect = c(low_trait, high_trait)), options = list() ))
ggplot(em45b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_negativeaffect)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))

#Detachment
#original model
#m46 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_detachment) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_detachment) +
#                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_detachment) + (1|id/run), trial_df)
#add random slopes
#m46 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_detachment) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_detachment) +
#                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_detachment) + (rt_lag_sc + reward_lag_rec|id) + (1|id/run), trial_df)
#add reward*rt_lag|ID
m46 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_detachment) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_detachment) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_detachment) + (rt_lag_sc*reward_lag_rec|id) + (1|id/run), trial_df)
#Add condition*WSLS 3-way and condition*trait 2-way interactions
#m46 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_detachment) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_detachment) +
#                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_detachment) + rewFunc*rt_lag_sc*reward_lag_rec +
#                rewFunc*scale(pid5_detachment) + (1|id/run), trial_df)
car::Anova(m46, '3')
summary(m46)
avg_trait <- mean(trial_df$pid5_detachment, na.rm=TRUE)
sd_trait <- sd(trial_df$pid5_detachment, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em46a <- as_tibble(emmeans::emtrends(m46, var = "rt_lag_sc", specs = c("reward_lag_rec", "pid5_detachment"), at = list(pid5_detachment = c(low_trait, high_trait))))
ggplot(em46a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_detachment)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) +  
    scale_y_reverse()
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em46b <- as_tibble(emmeans::emtrends(m46, var = "rt_vmax_lag_sc", specs = c("pid5_detachment", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), pid5_detachment = c(low_trait, high_trait)), options = list() ))
ggplot(em46b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_detachment)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))

#Antagonism 
#original model
#m47 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_antagonism) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_antagonism) +
#                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_antagonism) + (1|id/run), trial_df)
#add random slopes
#m47 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_antagonism) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_antagonism) +
#                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_antagonism) + (rt_lag_sc + reward_lag_rec|id) + (1|id/run), trial_df)
#add reward*rt_lag|ID
m47 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_antagonism) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_antagonism) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_antagonism) + (rt_lag_sc*reward_lag_rec|id) + (1|id/run), trial_df)
#Add condition*WSLS 3-way and condition*trait 2-way interactions
#m47 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_antagonism) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_antagonism) +
#                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_antagonism) + rewFunc*rt_lag_sc*reward_lag_rec +
#                rewFunc*scale(pid5_antagonism) + (1|id/run), trial_df)
car::Anova(m47, '3')
summary(m47)
avg_trait <- mean(trial_df$pid5_antagonism, na.rm=TRUE)
sd_trait <- sd(trial_df$pid5_antagonism, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em47a <- as_tibble(emmeans::emtrends(m47, var = "rt_lag_sc", specs = c("reward_lag_rec", "pid5_antagonism"), at = list(pid5_antagonism = c(low_trait, high_trait))))
ggplot(em47a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_antagonism)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) +  
    scale_y_reverse()
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em47b <- as_tibble(emmeans::emtrends(m47, var = "rt_vmax_lag_sc", specs = c("pid5_antagonism", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), pid5_antagonism = c(low_trait, high_trait)), options = list() ))
ggplot(em47b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_antagonism)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))

#Disinhibition
#original model
#m48 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_disinhibition) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_disinhibition) +
#                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_disinhibition) + (1|id/run), trial_df)
#add random slopes
#m48 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_disinhibition) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_disinhibition) +
#                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_disinhibition) + (rt_lag_sc + reward_lag_rec|id) + (1|id/run), trial_df)
#add reward*rt_lag|ID
m48 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_disinhibition) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_disinhibition) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_disinhibition) + (rt_lag_sc*reward_lag_rec|id) + (1|id/run), trial_df)
#Add condition*WSLS 3-way and condition*trait 2-way interactions
m48 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_disinhibition) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_disinhibition) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_disinhibition) + rewFunc*rt_lag_sc*reward_lag_rec +
                rewFunc*scale(pid5_disinhibition) + (1|id/run), trial_df)
car::Anova(m48, '3')
summary(m45)
avg_trait <- mean(trial_df$pid5_disinhibition, na.rm=TRUE)
sd_trait <- sd(trial_df$pid5_disinhibition, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em48a <- as_tibble(emmeans::emtrends(m48, var = "rt_lag_sc", specs = c("reward_lag_rec", "pid5_disinhibition"), at = list(pid5_disinhibition = c(low_trait, high_trait))))
ggplot(em48a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_disinhibition)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) +  
    scale_y_reverse()
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em48b <- as_tibble(emmeans::emtrends(m48, var = "rt_vmax_lag_sc", specs = c("pid5_disinhibition", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), pid5_disinhibition = c(low_trait, high_trait)), options = list() ))
ggplot(em48b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_disinhibition)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))

#Psychoticism 
#original model
#m49 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_psychoticism) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_psychoticism) +
#                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_psychoticism) + (1|id/run), trial_df)
#add random slopes
#m49 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_psychoticism) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_psychoticism) +
#                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_psychoticism) + (rt_lag_sc + reward_lag_rec|id) + (1|id/run), trial_df)
#add reward*rt_lag|ID
m49 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_psychoticism) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_psychoticism) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_psychoticism) + (rt_lag_sc*reward_lag_rec|id) + (1|id/run), trial_df)
#Add condition*WSLS 3-way and condition*trait 2-way interactions
#m49 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(pid5_psychoticism) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(pid5_psychoticism) +
#                rt_lag_sc*condition_trial_neg_inv_sc*scale(pid5_psychoticism) + rewFunc*rt_lag_sc*reward_lag_rec +
#                rewFunc*scale(pid5_psychoticism) + (1|id/run), trial_df)
car::Anova(m49, '3')
summary(m49)
avg_trait <- mean(trial_df$pid5_psychoticism, na.rm=TRUE)
sd_trait <- sd(trial_df$pid5_psychoticism, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em49a <- as_tibble(emmeans::emtrends(m49, var = "rt_lag_sc", specs = c("reward_lag_rec", "pid5_psychoticism"), at = list(pid5_psychoticism = c(low_trait, high_trait))))
ggplot(em49a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_psychoticism)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) +  
    scale_y_reverse()
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em49b <- as_tibble(emmeans::emtrends(m49, var = "rt_vmax_lag_sc", specs = c("pid5_psychoticism", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), pid5_psychoticism = c(low_trait, high_trait)), options = list() ))
ggplot(em49b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=pid5_psychoticism)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))



### COGNITIVE EFFECTS
#WTAR
m50 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(wtar) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(wtar) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(wtar) + (1|id/run), trial_df)
Anova(m50, '3')
summary(m50)
avg_trait <- mean(trial_df$wtar, na.rm=TRUE)
sd_trait <- sd(trial_df$wtar, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em50a <- as_tibble(emmeans::emtrends(m50, var = "rt_lag_sc", specs = c("reward_lag_rec", "wtar"), at = list(wtar = c(low_trait, high_trait))))
ggplot(em50a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=wtar)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) +  
    scale_y_reverse()
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em50b <- as_tibble(emmeans::emtrends(m50, var = "rt_vmax_lag_sc", specs = c("wtar", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), wtar = c(low_trait, high_trait)), options = list() ))
ggplot(em50b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=wtar)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))

#EXIT
m51 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(exit) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(exit) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(exit) + (1|id/run), trial_df)
Anova(m51, '3')
summary(m51)
avg_trait <- mean(trial_df$exit, na.rm=TRUE)
sd_trait <- sd(trial_df$exit, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em51a <- as_tibble(emmeans::emtrends(m51, var = "rt_lag_sc", specs = c("reward_lag_rec", "exit"), at = list(exit = c(low_trait, high_trait))))
ggplot(em51a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=exit)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) +  
    scale_y_reverse()
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em51b <- as_tibble(emmeans::emtrends(m51, var = "rt_vmax_lag_sc", specs = c("exit", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), exit = c(low_trait, high_trait)), options = list() ))
ggplot(em51b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=exit)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))


### Look at effect of Average Reward Rate
library(zoo)
library(car)
trial_df <- trial_df %>% group_by(id,run_number) %>% arrange(id, run_number, trial) %>% mutate(
    avr10 = scale(lag(rollmean(score_csv,10,na.pad = TRUE),10)), #AI added 12/5/22, average reward rate over previous 10 trials
    avr5 = scale(lag(rollmean(score_csv,5,na.pad = TRUE),5)), #AI added 12/5/22, average reward rate over previous 10 trials
    avr7 = scale(lag(rollmean(score_csv,7,na.pad = TRUE),7)), #AI added 12/5/22, average reward rate over previous 10 trials
    avr3 = scale(lag(rollmean(score_csv,3,na.pad = TRUE),3)), #AI added 12/5/22, average reward rate over previous 10 trials
) %>% ungroup()
#avr3
m60 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(avr3) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(avr3) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(avr3) + (1|id/run), trial_df)
car::Anova(m60, '3')
summary(m60)
avg_trait <- mean(trial_df$avr3, na.rm=TRUE)
sd_trait <- sd(trial_df$avr3, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em60a <- as_tibble(emmeans::emtrends(m60, var = "rt_lag_sc", specs = c("reward_lag_rec", "avr3"), at = list(avr3 = c(low_trait, high_trait))))
ggplot(em60a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=avr3)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) +  
    scale_y_reverse()
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em60b <- as_tibble(emmeans::emtrends(m60, var = "rt_vmax_lag_sc", specs = c("avr3", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), avr3 = c(low_trait, high_trait)), options = list() ))
ggplot(em60b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=avr3)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))
#avr5
m61 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(avr5) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(avr5) +
                rt_lag_sc*condition_trial_neg_inv_sc*scale(avr5) + (1|id/run), trial_df)
car::Anova(m61, '3')
summary(m61)
avg_trait <- mean(trial_df$avr5, na.rm=TRUE)
sd_trait <- sd(trial_df$avr5, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em61a <- as_tibble(emmeans::emtrends(m61, var = "rt_lag_sc", specs = c("reward_lag_rec", "avr5"), at = list(avr5 = c(low_trait, high_trait))))
ggplot(em61a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=avr5)) +
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) +  
    scale_y_reverse()
#Plot exploitation
early_trial <- -0.543 #all_data$condition_trial_neg_inv_sc[5], 5th trial (-0.5434399)
late_trial <- 0.479 #all_data$condition_trial_neg_inv_sc[40], 40th trial (+0.4787617)
em61b <- as_tibble(emmeans::emtrends(m61, var = "rt_vmax_lag_sc", specs = c("avr5", "condition_trial_neg_inv_sc"), at=list(condition_trial_neg_inv_sc = c(early_trial, late_trial), avr5 = c(low_trait, high_trait)), options = list() ))
ggplot(em61b, aes(x=condition_trial_neg_inv_sc, y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=avr5)) +
    geom_point(position = position_dodge(width = .6), size=2.5) +
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
    theme_bw(base_size=12) +ylab("Convergence on\n best RT (AU)")  +
    labs(shape = "Lethality") +
    theme(#axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + scale_x_discrete(name ="Trial", limits=c(-0.543, 0.479), labels=c("5","40"))

pdf("explore_RT_subj_byrewfunc_AV10.pdf", height = 15, width = 20)
ggplot(trial_df) + geom_point(aes(x=trial, y=rt_csv, color=rewFunc), size=1) + geom_point(aes(x=trial, y=avr10), size=0.75) + facet_wrap(~id)
dev.off()

anova_avr10 <- aov(avr10~rewFunc, data=trial_df, na.action=na.exclude)
summary(anova_avr10)
TukeyHSD(anova_avr10, na.action=na.exclude)


#See if RT variables differ between group
anova_two_way <- aov(rt_csv~groupLeth, data=trial_df)
summary(anova_two_way)
TukeyHSD(anova_two_way)

anova_two_way1 <- aov(rt_vmax~groupLeth, data=trial_df)
summary(anova_two_way1)
TukeyHSD(anova_two_way1)


#Sep 8, 2023. See if model parameters differ between groups.
#load demos
model_fits <- read.csv("/Users/angela/Documents/RESEARCH/ANALYSES/temporal_instrumental_agent/clock_task/vba_fmri/RESULTS/vba_out/explore/mfx/decay_factorize_selective_psequate_fixedparams_fmri/explore_decay_factorize_selective_psequate_fixedparams_fmri_mfx_sceptic_global_statistics.csv")
model_demos <- left_join(demos,model_fits,by="id")
summary(model_demos$Group)
summ_model <- model_demos %>% group_by(Group) %>%
    summarise(group_n = n(),
              alpha_mean = mean(alpha, na.rm=TRUE),
              alpha_sd = sd(alpha, na.rm=TRUE),
              beta_mean = mean(beta, na.rm=TRUE),
              beta_sd = sd(beta, na.rm=TRUE),
              gamma_mean = mean(gamma, na.rm=TRUE),
              gamma_sd = sd(gamma,na.rm=TRUE),
              r2_mean = mean(as.numeric(R2),na.rm=TRUE),
              r2_sd = sd(as.numeric(R2),na.rm=TRUE))
print(summ_model)

#Run statistical tests 
alpha_aov <- aov(alpha ~ Group, data = model_demos)
summary(alpha_aov)
beta_aov <- aov(beta ~ Group, data = model_demos)
summary(beta_aov)
gamma_aov <- aov(gamma ~ Group, data = model_demos)
summary(gamma_aov)
r2_aov <- aov(as.numeric(R2) ~ Group, data = model_demos)
summary(r2_aov)

