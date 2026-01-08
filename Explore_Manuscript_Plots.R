# Plots for Explore paper
# Dec 2025
# Angela Ianni

rm(list=ls())

library(stringr)
library(pracma)
library(wesanderson)
library(tidyverse)
library(data.table) 
library(formula.tools)
library(RColorBrewer)

########## User Settings ##########
subjs_to_run = "attempters_only" #options are all, controls_only, patients_only, or attempters_only
trial_data = ("/Users/angela/Documents/Work/Explore/Medusa_analysis/trial_df_03282023.rds")
subject_demographics = ("/Users/angela/Documents/Work/Explore/Medusa_analysis/explore_complete_demographics.rds") #Updated Nov 2025 with fixed impulsivity measures
excludes = c(207224, 210374, 210701, 213708, 216806, 216845, 219392, 220024, 221698, 440149, 440311, 440336) #imaging QC excludes (excluded in original anlaysis); excludes = c(207224, 210374, 210701, 213708, 216806, 216845, 219392, 220024, 221698, 440149, 440311, 440336, 211253, 431100, 440243, 440254, 440263, 440392, 221292, 221637, 220913, 431224) #imaging QC excludes, imaging questional QCs, behavioral flatliners, 2 subjs with really long runs
align="clock" #options: clock, feedback

########## Load trial and demographic data ##########
trial_df <- readRDS(trial_data)
trial_df$id <- as.integer(sub("_1", "", trial_df$id))
load('/Users/angela/Documents/Work/Explore/Medusa_analysis/combined_rs_betas.Rda')
# Recode omissions as -0.5 and rewards as 0.5 and code the condition trial (more similar to bsocial run_trial)
trial_df <- trial_df %>% 
  group_by(id,run_number) %>% 
  arrange(id, run_number, trial) %>% 
  mutate(reward_lag_rec = if_else(reward_lag=="omission", -0.5, 0.5),
         condition_trial = run_trial-floor(run_trial/40.5)*40,                                                                                          
         condition_trial_neg_inv = -1000 / condition_trial,
         condition_trial_neg_inv_sc = as.vector(scale(condition_trial_neg_inv)),
  ) %>% ungroup()
demos <- readRDS(subject_demographics)

#Fix mis-named IDs for Explore
trial_df$id <- gsub(x = trial_df$id, pattern='881224',replacement='431224')
trial_df$id <- gsub(x = trial_df$id, pattern='881230',replacement='431230')


# convert id to character
trial_df$id <- as.numeric(trial_df$id)

#Filter out excludes
trial_df <- trial_df %>% filter(!(id %in% excludes)) 
demos <- demos %>% filter(!(id %in% excludes)) 

#Merge demos and trial_df
trial_df <- trial_df %>% left_join(demos, by="id") 
#Merge WSLS behavior and trial_df
wsls_df <- combined_rs_betas %>% dplyr::select(id,WSLS_rs,LS_rs)
wsls_df$id <- as.integer(wsls_df$id)
demos_wsls <- demos %>% left_join(wsls_df, by="id")

demos_wsls$groupLeth <- factor(demos_wsls$groupLeth, levels = c("HL_Attempters", "LL_Attempters", "Ideators", "Depressed", "Controls"))


#ggplot(demos_wsls, aes(x = groupLeth, y = LS_rs)) +
#  +     stat_summary(
#    +         fun = mean,
#    +         fun.data = mean_se,
#    +         geom = "pointrange",
#    +         color = "black"
#    +     ) +
#  +     theme_classic() +
#  +     labs(
#    +         x = "Lethality Group",
#    +         y = "LS_rs (mean ± SE)"
#    +     )

#my_cols <- c(
#  "HL_Attempters" = "",
#  "LL_Attempters" = "",
#  "Ideators" = "",
#  "Depressed" = "",
#  "Controls" = ""
#)


### Plot WSLS by group
# Define the ColorBrewer palette
brewer_palette_tmp <- brewer.pal(n = 5, name = "Set1")
brewer_palette <- brewer_palette_tmp
brewer_palette[3] <- brewer_palette_tmp[5]
brewer_palette[4] <- brewer_palette_tmp[3]
brewer_palette[5] <- brewer_palette_tmp[4]

ggplot(demos_wsls, aes(x = groupLeth, y = WSLS_rs, color = groupLeth)) +
       stat_summary(
             fun = mean,
             fun.data = mean_se,
             geom = "pointrange",
             size = 1.5,
             linewidth = 1.5
         ) +
       theme_classic() +
  theme(axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=16, lineheight = 1.1, vjust = 0.5)
  ) +
  scale_x_discrete(labels = NULL) +
  scale_color_manual(values = brewer_palette) +
       labs(
             x = "Lethality Group",
             y = "WSLS random slopes\n(mean ± SE)"
         )

#ANOVA for differences in WSLS random slopes
aov_wsls_groupLeth <- aov(WSLS_rs ~ groupLeth, data=demos_wsls)
aov_wsls_Group <- aov(WSLS_rs ~ Group, data=demos_wsls)
summary(aov_wsls_groupLeth)
TukeyHSD(aov_wsls_groupLeth)
summary(aov_wsls_Group)
TukeyHSD(aov_wsls_Group)


trial_df$groupLeth <- factor(trial_df$groupLeth, levels = c("HL_Attempters", "LL_Attempters", "Ideators", "Depressed", "Controls"))

#### Behavioral Models
#orig model
#m1 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*groupLeth + rt_vmax_lag_sc*condition_trial_neg_inv*groupLeth + rt_lag_sc*condition_trial_neg_inv*groupLeth + (1|id/run), trial_df)
m1 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*groupLeth + rt_vmax_lag_sc*trial_neg_inv*groupLeth + rt_lag_sc*trial_neg_inv*groupLeth + (1|id/run), trial_df)
car::Anova(m1,'3')
#summary(m1)
m1a <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*Group + rt_vmax_lag_sc*trial_neg_inv*Group + rt_lag_sc*trial_neg_inv*Group + (1|id/run), trial_df)
car::Anova(m1a,'3')
#plot GroupLeth
em1 <- as_tibble(emmeans::emtrends(m1, var = "rt_lag_sc", specs = c("reward_lag_rec", "groupLeth")))
em1$reward_lag_rec <- factor(em1$reward_lag_rec, levels = c("-0.5", "0.5"), labels = c("Omission", "Reward"))
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
ggplot(em1, aes(x = reward_lag_rec, y = rt_lag_sc.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = groupLeth)) +
  geom_point(position = position_dodge(width = .6), size = 2.5) + 
  geom_errorbar(position = position_dodge(width = 0.6), width = 0.4, linewidth = 0.9) + 
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

#Pairwise comparisons of groupLeth (Added 10/21/2024 - don't think this is right)
#pair_m1.emm <- emmeans(m1, pairwise ~ reward_lag_rec*groupLeth)
#pairs(pair_m1.emm, adjust="fdr")

#Plot Group (not separating lethality groups)
em1a <- as_tibble(emmeans::emtrends(m1a, var = "rt_lag_sc", specs = c("reward_lag_rec", "Group")))
em1a$reward_lag_rec <- factor(em1a$reward_lag_rec, levels = c("-0.5", "0.5"), labels = c("Omission", "Reward"))
em1a$Group <- factor(em1a$Group, level = c("Attempters", "Ideators", "Depressed", "Controls"))
ggplot(em1a, aes(x = reward_lag_rec, y = rt_lag_sc.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = Group)) +
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
                     labels = c("Attempters", "Ideators", "Depressed", "Controls"))





#### Load in the GLM betas
fmri_dir <- "~/Documents/Work/Explore/fMRI/Dec2022/extracted_values/12Dec2022/"
setwd(fmri_dir)

# read in V_max betas (or zstat or tstat or cope)
if (align=="clock") {
  setwd(file.path(fmri_dir, "L1m-v_max_wi")) #clock-aligned
  vmax <- read_csv("transformed_Schaefer_244_final_3.125mm_cope_l2.csv.gz") %>% filter(l2_model == "intercept" & l1_cope_name == "EV_v_max_wi" & !(id %in% excludes))
} else if (align=="feedback") {
  setwd(file.path(fmri_dir, "L1m-v_max_wi_lead")) #Jan 2023, feedback aligned
  vmax <- read_csv("transformed_Schaefer_244_final_3.125mm_cope_l2.csv.gz") %>% filter(l2_model == "intercept" & l1_cope_name == "EV_v_max_wi_lead" & !(id %in% excludes))
}
setwd(fmri_dir)
labels <- read_csv("region_labels_244.csv") %>% mutate(roi_num = as.numeric(roi_num)) %>% inner_join(read_csv("region_lookup_244.csv"), by = "roi_num")

vmax <- vmax %>% mutate(roi_num = mask_value) %>% dplyr::select(-x, -y, -z) %>% inner_join(labels, by = "roi_num") %>% inner_join(demos, by = "id")

# label missing networks as subcortical
vmax$network[is.na(vmax$network)] <- "subcortical"

### April 15, 2024 - check for betas association with behavior
vPFC <- vmax %>% 
  filter(roi_num %in% c(84,86,88,89,192,194,65,66,67,170,171,55,56,159,160,161,191)) %>% dplyr::select(c(network, subregion, value, hemi, id)) %>%  mutate(region_network_label = paste(network, subregion, hemi, sep = "_")) 
#sub <- vmax %>% #only includes amygdala, hippocampus, DMN, Limbic, and control networks of striatum
#  filter(roi_num %in% c(201,202,203,204,215,216,217,220,221,222,223,224,225,226,227,228,229,230)) %>% 
#  select(c(network, subregion, value, hemi, id)) %>%  mutate(region_network_label = paste(network, subregion, hemi, sep = "_")) 
sub <- vmax %>% # includes amygdala, hippocampus, and all networks of striatum (adjusted 4/17/2024)
  filter(roi_num %in% c(201,202,203,204,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230)) %>% 
  dplyr::select(c(network, subregion, value, hemi, id)) %>%  mutate(region_network_label = paste(network, subregion, hemi, sep = "_"))
vPFC_DMN <- vPFC %>% filter(network=="Default")
vPFC_LIM <- vPFC %>% filter(network=="Limbic")
vPFC_CON <- vPFC %>% filter(network=="Cont")
sub_DMN <- sub %>% filter(network=="Default")
sub_LIM <- sub %>% filter(network=="Limbic")
sub_CON <- sub %>% filter(network=="Cont")
sub_SalVentAttn <- sub %>% filter(network=="SalVentAttn")
sub_SomMot <- sub %>% filter(network=="SomMot")
sub_HC <- sub %>% filter(str_detect(subregion, regex("Hippocampus", ignore_case=TRUE)))
sub_amyg <- sub %>% filter(str_detect(subregion, regex("BLA|CMN", ignore_case=TRUE)))

#Find average beta within each network for each subject
vPFC_DMN_mean <- vPFC_DMN %>% group_by(id) %>%
  dplyr::summarize(vPFC_DMN_mean=mean(value,na.rm=TRUE)) %>% ungroup()
vPFC_LIM_mean <- vPFC_LIM %>% group_by(id) %>%
  dplyr::summarize(vPFC_LIM_mean=mean(value,na.rm=TRUE)) %>% ungroup()
vPFC_CON_mean <- vPFC_CON %>% group_by(id) %>%
  dplyr::summarize(vPFC_CON_mean=mean(value,na.rm=TRUE)) %>% ungroup()
sub_DMN_mean <- sub_DMN %>% group_by(id) %>%
  dplyr::summarize(sub_DMN_mean=mean(value,na.rm=TRUE)) %>% ungroup()
sub_LIM_mean <- sub_LIM %>% group_by(id) %>%
  dplyr::summarize(sub_LIM_mean=mean(value,na.rm=TRUE)) %>% ungroup()
sub_CON_mean <- sub_CON %>% group_by(id) %>%
  dplyr::summarize(sub_CON_mean=mean(value,na.rm=TRUE)) %>% ungroup()
sub_SalVentAttn_mean <- sub_SalVentAttn %>% group_by(id) %>%
  dplyr::summarize(sub_SalVentAttn_mean=mean(value,na.rm=TRUE)) %>% ungroup()
sub_SomMot_mean <- sub_SomMot %>% group_by(id) %>%
  dplyr::summarize(sub_SomMot_mean=mean(value,na.rm=TRUE)) %>% ungroup()
sub_HC_mean <- sub_HC %>% group_by(id) %>%
  dplyr::summarize(sub_HC_mean=mean(value,na.rm=TRUE)) %>% ungroup()
sub_amyg_mean <- sub_amyg %>% group_by(id) %>%
  dplyr::summarize(sub_amyg_mean=mean(value,na.rm=TRUE)) %>% ungroup()

#trial_df <- trial_df %>%
#  left_join(demos, by="id") 
#trial_df <- trial_df %>% rename(groupLeth=group_leth,sex=gender)
trial_df <- trial_df %>%
  left_join(vPFC_DMN_mean, by="id") 
trial_df <- trial_df %>%
  left_join(vPFC_LIM_mean, by="id") 
trial_df <- trial_df %>%
  left_join(vPFC_CON_mean, by="id") 
trial_df <- trial_df %>%
  left_join(sub_DMN_mean, by="id")
trial_df <- trial_df %>%
  left_join(sub_LIM_mean, by="id") 
trial_df <- trial_df %>%
  left_join(sub_CON_mean, by="id") 
trial_df <- trial_df %>%
  left_join(sub_SalVentAttn_mean, by="id") 
trial_df <- trial_df %>%
  left_join(sub_SomMot_mean, by="id") 
trial_df <- trial_df %>%
  left_join(sub_HC_mean, by="id")
trial_df <- trial_df %>%
  left_join(sub_amyg_mean, by="id")
trial_df$groupLeth <- factor(trial_df$groupLeth, level = c("HL_Attempters", "LL_Attempters", "Ideators", "Depressed", "Controls"))

LS_trial_df <- trial_df %>% filter(reward_lag_rec==-0.5)

#add WSLS_rs and LS_rs to trial_df
trial_df <- trial_df %>% left_join(wsls_df, by="id")
demos_wsls <- demos_wsls %>% left_join(vPFC_DMN_mean, by = "id")





trial_df_att <- trial_df %>% filter(Group == "Attempters")
demos_wsls_att <- demos_wsls %>% filter(Group == "Attempters")
###Max intent
m5 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(max_intent) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(max_intent) +
             rt_lag_sc*condition_trial_neg_inv_sc*scale(max_intent) + (1|id/run), trial_df_att)
car::Anova(m5, type=3)
summary(m5)
avg_trait <- mean(trial_df_att$max_intent, na.rm=TRUE)
sd_trait <- sd(trial_df_att$max_intent, na.rm=TRUE)
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
###Max Lethality
m6 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(max_lethality) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(max_lethality) +
             rt_lag_sc*condition_trial_neg_inv_sc*scale(max_lethality) + (1|id/run), trial_df_att)
car::Anova(m6, type=3)
summary(m6)
avg_trait <- mean(trial_df_att$max_lethality, na.rm=TRUE)
sd_trait <- sd(trial_df_att$max_lethality, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em6a <- as_tibble(emmeans::emtrends(m6, pbkrtest.limit = 7767, var = "rt_lag_sc", specs = c("reward_lag_rec", "max_lethality"), at = list(max_lethality = c(low_trait, high_trait))))
ggplot(em6a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=lower.CL, ymax=upper.CL, color=max_lethality)) +
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) +  
  scale_y_reverse()
#Max lethality, updated format:
avg_trait <- mean(trial_df_att$max_lethality, na.rm = TRUE)
sd_trait  <- sd(trial_df_att$max_lethality, na.rm = TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait<- avg_trait + 2*sd_trait

em6a <- as_tibble(emmeans::emtrends(m6, pbkrtest.limit = 7767, var = "rt_lag_sc",
                                    specs = c("reward_lag_rec","max_lethality"),
                                    at = list(max_lethality = c(low_trait, high_trait))))

em6a$reward_lag_rec <- factor(em6a$reward_lag_rec, levels = c("-0.5","0.5"),
                              labels = c("Omission","Reward"))
em6a$max_lethality  <- factor(em6a$max_lethality,
                              levels = c(low_trait, high_trait),
                              labels = c("Low","High"))

pd <- position_dodge(width = 0.3)

ggplot(em6a, aes(x = reward_lag_rec, y = rt_lag_sc.trend,
                 ymin = lower.CL, ymax = upper.CL,
                 color = max_lethality, group = max_lethality)) +
  geom_errorbar(position = pd, width = 0.3, linewidth = 1.5) +
  geom_point(position = pd, size = 3) +
  theme_bw(base_size = 12) +
  labs(x = NULL, y = "RT swings\n Small <---------> Large", color = "Max lethality") +
  scale_color_manual(values = c("Low" = "#4DAF4A", "High" = "#984EA3")) +
  scale_x_discrete(drop = FALSE) +
  scale_y_reverse() +
  theme(axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12, color = "grey10"),
        panel.grid.major.x = element_blank(),
        legend.text  = element_text(size = 12),
        legend.title = element_text(size = 14))
###Max Ideation
m7 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(max_ideation) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(max_ideation) +
             rt_lag_sc*condition_trial_neg_inv_sc*scale(max_ideation) + (1|id/run), trial_df_att)
car::Anova(m7, type=3)
summary(m7)
avg_trait <- mean(trial_df_att$max_ideation, na.rm=TRUE)
sd_trait <- sd(trial_df_att$max_ideation, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em7a <- as_tibble(emmeans::emtrends(m7, pbkrtest.limit = 7767, var = "rt_lag_sc", specs = c("reward_lag_rec", "max_ideation"), at = list(max_ideation = c(low_trait, high_trait))))
ggplot(em7a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=lower.CL, ymax=upper.CL, color=max_ideation)) +
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) +  
  scale_y_reverse()
###Age at first attempt
m8 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(age_at_first_attempt) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(age_at_first_attempt) +
             rt_lag_sc*condition_trial_neg_inv_sc*scale(age_at_first_attempt) + (1|id/run), trial_df_att)
car::Anova(m8, type=3)
summary(m8)
avg_trait <- mean(trial_df_att$age_at_first_attempt, na.rm=TRUE)
sd_trait <- sd(trial_df_att$age_at_first_attempt, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em8a <- as_tibble(emmeans::emtrends(m8, pbkrtest.limit = 7767, var = "rt_lag_sc", specs = c("reward_lag_rec", "age_at_first_attempt"), at = list(age_at_first_attempt = c(low_trait, high_trait))))
ggplot(em8a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=lower.CL, ymax=upper.CL, color=age_at_first_attempt)) +
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) +  
  scale_y_reverse()
###vPFC DMN Vmax
m9 <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec*scale(vPFC_DMN_mean) + rt_vmax_lag_sc*condition_trial_neg_inv_sc*scale(vPFC_DMN_mean) +
             rt_lag_sc*condition_trial_neg_inv_sc*scale(vPFC_DMN_mean) + (1|id/run), trial_df)
car::Anova(m9, type=3)
summary(m9)
avg_trait <- mean(trial_df$vPFC_DMN_mean, na.rm=TRUE)
sd_trait <- sd(trial_df$vPFC_DMN_mean, na.rm=TRUE)
low_trait <- avg_trait - 2*sd_trait
high_trait <- avg_trait + 2*sd_trait
em9a <- as_tibble(emmeans::emtrends(m9, pbkrtest.limit = 7767, var = "rt_lag_sc", specs = c("reward_lag_rec", "vPFC_DMN_mean"), at = list(vPFC_DMN_mean = c(low_trait, high_trait))))
ggplot(em9a, aes(x=reward_lag_rec, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=vPFC_DMN_mean)) +
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  theme_bw(base_size=12) +ylab("RT swings (AU)\n Small <---------> Large")  + 
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) +  
  scale_y_reverse()




# Plot the WSLS effects on a single figure
em_leth <- em6a %>% mutate (
  var = case_when(
  #max_lethality_SD = case_when(
    max_lethality < "0" ~ "Low",
    max_lethality > "7" ~ "High")
) %>% rename(
  max_lethality.rt_lag_sc.trend = rt_lag_sc.trend,
  max_lethality.SE = SE,
  max_lethality.df = df,
  max_lethality.lower.CL = lower.CL,
  max_lethality.upper.CL = upper.CL
) %>% dplyr::select(!max_lethality)

em_intent <- em5a %>% mutate (
  var = case_when(
  #max_intent_SD = case_when(
    max_intent < "10" ~ "Low",
    max_intent > "20" ~ "High")
) %>% rename(
  max_intent.rt_lag_sc.trend = rt_lag_sc.trend,
  max_intent.SE = SE,
  max_intent.df = df,
  max_intent.lower.CL = asymp.LCL,
  max_intent.upper.CL = asymp.UCL
) %>% dplyr::select(!max_intent)

em_agefirst <- em8a %>% mutate (
  #age_at_first_attempt_SD = case_when(
  var = case_when(
    age_at_first_attempt < "0" ~ "Low",
    age_at_first_attempt > "70" ~ "High")
) %>% rename(
  age_at_first_attempt.rt_lag_sc.trend = rt_lag_sc.trend,
  age_at_first_attempt.SE = SE,
  age_at_first_attempt.df = df,
  age_at_first_attempt.lower.CL = lower.CL,
  age_at_first_attempt.upper.CL = upper.CL
) %>% dplyr::select(!age_at_first_attempt)

em_ideation <- em7a %>% mutate (
  #max_ideation_SD = case_when(
  var = case_when(
    max_ideation < "20" ~ "Low",
    max_ideation > "40" ~ "High")
) %>% rename(
  max_ideation.rt_lag_sc.trend = rt_lag_sc.trend,
  max_ideation.SE = SE,
  max_ideation.df = df,
  max_ideation.lower.CL = lower.CL,
  max_ideation.upper.CL = upper.CL
) %>% dplyr::select(!max_ideation)

em_vPFCDMN <- em9a %>% mutate (
  var = case_when(
  #vPFC_DMN_mean_SD = case_when(
    vPFC_DMN_mean < "0" ~ "Low",
    vPFC_DMN_mean > "0" ~ "High")
) %>% rename(
  vPFC_DMN_mean.rt_lag_sc.trend = rt_lag_sc.trend,
  vPFC_DMN_mean.SE = SE,
  vPFC_DMN_mean.df = df,
  vPFC_DMN_mean.lower.CL = asymp.LCL,
  vPFC_DMN_mean.upper.CL = asymp.UCL
) %>% dplyr::select(!vPFC_DMN_mean)

list_of_dfs <- list(em_leth, em_intent, em_ideation, em_agefirst, em_vPFCDMN)
common_cols <- c("reward_lag_rec", "var")
#combined_ems <- merge(em_leth, em_intent, em_ideation, em_agefirst, em_vPFCDMN, by = c("reward_lag_rec", "var"))

combined_df <- reduce(list_of_dfs, inner_join, by = common_cols)

#### Plot them all on a single plot
# 1. Reshape combined_df to long tidy format
plot_df <- combined_df %>% 
  pivot_longer(
    # keep x-variable (and any grouping var you want)
    cols = -c(reward_lag_rec, var),
    names_to = c("measure", "stat"),
    names_pattern = "^(.*)\\.(rt_lag_sc\\.trend|lower\\.CL|upper\\.CL)$"
  ) %>%
  # 2. Spread stat columns back out
  pivot_wider(
    names_from  = stat,
    values_from = value
  )

# plot_df now has:
# reward_lag_rec, var, measure,
# `rt_lag_sc.trend`, `lower.CL`, `upper.CL`
# 3. Plot all measures using your original style
ggplot(
  plot_df,
  aes(
    x    = reward_lag_rec,
    y    = `rt_lag_sc.trend`,
    ymin = `lower.CL`,
    ymax = `upper.CL`,
    color = measure          # or color = var, if that makes more sense
  )
) +
  geom_point(position = position_dodge(width = 0.6), size = 2.5) +
  geom_errorbar(position = position_dodge(width = 0.6),
                width = 0.4, size = 0.9) +
  theme_bw(base_size = 12) +
  ylab("RT swings (AU)\n Small <---------> Large") +
  theme(
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text = element_text(size = 8.5, color = "grey10")
  ) #+
  #scale_y_reverse() #+
  #facet_wrap(~ measure)   # remove this if you prefer one combined panel




#### See if vPFC DMN Vmax betas mediate the lethality * WSLS effect
# Load packages (run this at the start of each R session)
library(lavaan)
library(semPlot)

demos_wsls_att$max_lethality <- scale(demos_wsls_att$max_lethality)
demos_wsls_att$age_at_first_attempt <- scale(demos_wsls_att$age_at_first_attempt)
demos_wsls_att$max_intent <- scale(demos_wsls_att$max_intent)
demos_wsls_att$vPFC_DMN_mean <- scale(demos_wsls_att$vPFC_DMN_mean)
demos_wsls_att$WSLS_rs <- scale(demos_wsls_att$WSLS_rs)

#cor(demos_wsls_att[,c("age_at_first_attempt","vPFC_DMN_mean","WSLS_rs")], use="pairwise") #correlations are 0.22 w/ vPFC DMN mean, -0.32 w/ WSLS_rs
#cor(demos_wsls_att[,c("max_lethality","vPFC_DMN_mean","WSLS_rs")], use="pairwise") #correlations are -0.07 w/ vPFC DMN mean, 0.056 w/ WSLS_rs
#cor(demos_wsls_att[,c("max_intent","vPFC_DMN_mean","WSLS_rs")], use="pairwise") #correlations are -0.12 w/ vPFC DMN mean, -0.086 w/ WSLS_rs
cor(demos_wsls_att[,c("age_at_first_attempt","max_lethality","max_intent","vPFC_DMN_mean","WSLS_rs")], use="pairwise") #correlations are -0.12 w/ vPFC DMN mean, -0.086 w/ WSLS_rs

#Regression with all of these
model_vpfc <- lm(
  vPFC_DMN_mean ~ age_at_first_attempt + max_lethality + max_intent,
  data = demos_wsls_att
)
summary(model_vpfc)

model_wsls <- lm(
  WSLS_rs ~ age_at_first_attempt + max_lethality + max_intent,
  data = demos_wsls_att
)
summary(model_wsls)

model_vpfc_wsls_lethality <- lm(
  #max_lethality ~ vPFC_DMN_mean + WSLS_rs,
  max_lethality ~ vPFC_DMN_mean*WSLS_rs,
  data = demos_wsls_att
)
summary(model_vpfc_wsls_lethality)

model_vpfc_wsls_intent <- lm(
  #max_intent ~ vPFC_DMN_mean + WSLS_rs,
  max_intent ~ vPFC_DMN_mean*WSLS_rs,
  data = demos_wsls_att
)
summary(model_vpfc_wsls_intent)

model_vpfc_wsls_attemptage <- lm(
  #age_at_first_attempt ~ vPFC_DMN_mean + WSLS_rs,
  age_at_first_attempt ~ vPFC_DMN_mean*WSLS_rs,
  data = demos_wsls_att
)
summary(model_vpfc_wsls_attemptage)

#mediation model: max_lethality -> vPFC_DMN_mean -> WSLS_rs
med_model <- '
  #path a
  vPFC_DMN_mean ~ a*age_at_first_attempt
  
  #path b
  #WSLS_rs ~ b*vPFC_DMN_mean
  
  #direct effect
  WSLS_rs ~ b*vPFC_DMN_mean + c_prime*age_at_first_attempt
  
  #indirect effect
  ind := a*b
  
  #total effect (direct + indirect)
  total := c_prime + (a*b)
'


#### Multilevel mediation
trial_df_mediation <- trial_df
trial_df_mediation <- trial_df_mediation %>%
  mutate(interaction = rt_lag_sc * reward_lag_rec)
subject_level <- trial_df_mediation %>%
  group_by(id) %>%
  summarise(
    vPFC_DMN_mean = mean(vPFC_DMN_mean, na.rm = TRUE),
    interaction_mean = mean(interaction, na.rm = TRUE)#,
    #by = "id",
    #.groups = "drop"
  )

med_model <- lm(vPFC_DMN_mean ~ interaction_mean, data=subject_level)
subject_level_interaction <- dplyr::select(subject_level, id, interaction_mean)
trial_df_mediation <- trial_df_mediation %>% left_join(subject_level_interaction, by="id")
outcome_model <- lmer(rt_csv ~ rt_lag_sc*reward_lag_rec + interaction_mean + vPFC_DMN_mean + (1+rt_lag_sc | id), data=trial_df_mediation)
summary(outcome_model)

library(mediation)

med_out <- mediate(
  model.m = med_model,
  model.y = outcome_model,
  treat = "interaction_mean",
  mediator = "vPFC_DMN_mean",
  boot = FALSE, #default is quasi-Bayesian Monte Carlo simulation to get CIs
  sims = 5000
)
summary(med_out) #n.s.




#### Plot suicide variables
demos_att <- demos %>% dplyr::filter(Group == "Attempters")
ggplot(demos_att, aes(x = age_at_first_attempt, y = max_lethality)) +
       geom_point() +
       geom_smooth(method = "lm", se = TRUE) +
       labs(
           x = "Age at First Attempt",
             y = "Max Lethality",
             title = "Relationship Between Age at First Attempt and Max Lethality"
         ) +
      theme_minimal()

ggplot(demos_att, aes(x = age_at_first_attempt, y = max_intent)) +
       geom_point() +
       geom_smooth(method = "lm", se = TRUE) +
       labs(
             x = "Age at First Attempt",
             y = "Max Intent",
             title = "Relationship Between Age at First Attempt and Max Intent"
       ) +
       theme_minimal()

ggplot(demos_att, aes(x = max_lethality, y = max_intent)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "Max Lethality",
    y = "Max Intent",
    title = "Relationship Between Max Lethality and Max Intent"
  ) +
  theme_minimal()

