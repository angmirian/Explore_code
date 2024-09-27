# aux file containing commonly used functions


#### DATA MANAGEMENT ####
#recode pid5 to start with "pid5", and compute pid5 subscales
AI_recodePID5 <- function(data) {
  data <- data %>%
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
  data<-data %>% dplyr::mutate(pid5_negativeaffect = (pid5_anxiousness+ pid5_emotional_lability+ pid5_separation_insecurity)/3)
  data<-data %>% dplyr::mutate(pid5_detachment = (pid5_withdrawal+ pid5_anhedonia+ pid5_intimacy_avoidance)/3)
  data<-data %>% dplyr::mutate(pid5_antagonism = (pid5_manipulativeness+ pid5_deceitfulness+ pid5_grandiosity)/3)
  data<-data %>% dplyr::mutate(pid5_disinhibition = (pid5_irresponsibility+ pid5_impulsivity+ pid5_distractibility)/3)
  data<-data %>% dplyr::mutate(pid5_psychoticism = (pid5_unusual_beliefs_experiences+ pid5_eccentricity+ pid5_perceptual_dysregulation)/3)
}

#code simple groups (patients vs. controls)
AI_simplegroups <- function(data) {
  data$group_simple <- data$groupLeth
  data$group_simple[str_detect(data$group_simple, regex("Attempters", ignore_case=TRUE))] <- "Depressed"
  data$group_simple[str_detect(data$group_simple, regex("LL_Attempters", ignore_case=TRUE))] <- "Depressed"
  data$group_simple[str_detect(data$group_simple, regex("Ideators", ignore_case=TRUE))] <- "Depressed"
  data$group_simple <- recode_factor(data$group_simple, Controls = "HC", Depressed = "Depressed")
  data$group_simple <- droplevels(data$group_simple)
  data$group_simple <- factor(data$group_simple, level = c("HC","Depressed"))
}



#### Plotting functions ####
AI_histogram <- function(data, x_var, group_var, alpha_var, colors, bin_var) {
  data|>
    ggplot(aes(x=.data[[x_var]], color=.data[[group_var]], fill=.data[[group_var]])) +
    geom_histogram(alpha=alpha_var, position='identity', show.legend=FALSE, bins=bin_var) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors)
}

AI_boxplot <- function(data, y_var, group_var, alpha_var, colors) {
  data|>
    ggplot(aes(y=.data[[y_var]], x=.data[[group_var]], color=.data[[group_var]], fill=.data[[group_var]])) +
    geom_boxplot(alpha=alpha_var, show.legend=FALSE) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors)
}

