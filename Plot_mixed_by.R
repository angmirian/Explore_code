#AI - Plot_mixed_byI.R; uses ddf generated from Build_mixed_by_modular.R
#Last updated May 1, 2024

rm(list=ls())

library(wesanderson)
library(pracma)
library(ggplot2)
library(tidyverse)
library(viridis)
library(grid)
library(viridis) 

########## User Settings ##########
align_to = "clock" #options are clock or response 
structure = "network" 
toplot <- "estimate" #options: are estimate and statistic
region_to_run = "cortical" #options are cortical (vmPFC) and subcortical
subjs_to_run = "all" #options are all or controls_only
flip_stat = 1 #to flip statistic for plots (i.e. if model set up to look at attempters vs. depressed but you want to flip to see depressed vs. attempters)
if (align_to=="clock") {
  exclude_prefeedback = 0 #set to 1 if you want to exclude the pre-feedback results from plots
} else if (align_to=="response") {exclude_prefeedback = 1}
select_networks=TRUE #can only set to TRUE for subcortical; if true only do Cont, Default, Limbic, hipp, and amygdala

subcortical_cache_dir = '/Users/angela/Documents/Research/Explore/Medusa_analysis/subcortical' 
cortical_cache_dir = '/Users/angela/Documents/Research/Explore/Medusa_analysis/cortical'

if (region_to_run=="cortical") {data_subfolder <- "Patients_only"
} else {data_subfolder <- "Selected_subcortical_regions/Patients_only" }



########## Load data ##########
#Define where data is located and where plots will be saved
if(region_to_run == "subcortical") {
  if (structure == "network") {dir <- paste0(subcortical_cache_dir,'/Network_CombinedHemis/',data_subfolder)}
} else if(region_to_run == "cortical") {
  if (structure == "network") {dir <- paste0(cortical_cache_dir,'/Network_CombinedHemis/',data_subfolder)}
}
setwd(dir)
#Read list of model files
if (align_to=="clock") {
  #if (subjs_to_run=="all") {
    data_files <- read.table('clock_data_list.txt')
    #data_files <- read.table('clock_data_list.txt')
  #} else if (subjs_to_run=="controls_only") {
  #  data_files <- read.table('clock_data_list_CONTROLS_ONLY.txt')
  #}
} else if (align_to=="response") {
  #if (subjs_to_run=="all") {
      data_files <- read.table('response_data_list.txt')
      #data_files <- read.table('response_data_list.txt')
  #} else if (subjs_to_run=="controls_only") {
  #  data_files <- read.table('response_data_list_CONTROLS_ONLY.txt')
  #}
}



########## Loop through models ##########
for (i in 1:length(data_files$V1)) {
  model <- data_files$V1[i]
  print(model)
  model_name <- str_split(model,'.Rdata')[[1]][1]
  #Load the Rdata file with mixed_by results
  load(model)
  #Define epoch label
  if (align_to == "response"){epoch_label <- 'response'
  } else if(align_to == "clock"){epoch_label <- 'clock'}
  #Select what you want to plot
  if (toplot == "statistic") {ylabel <- 't-statistic'
  } else if(toplot == "estimate") {ylabel <- 'Coefficient (AU)'}
  
  ########## Get data ready to plot ##########
  # ddf is the output structure; reml and ml are two different convergence methods (use reml - contains AIC and BIC andn log likelihood)
  # ddf$coeff_df_reml - has all fixed and random effects; value is the estimate and estimated error is std.error
  #consolidate the ddf to just the coef_df_reml
  ddf <- ddf$coef_df_reml
  if (structure == "network") {ddf$network <- as.factor(ddf$network)}
  if (exclude_prefeedback) {
    ddf<- ddf %>% filter(evt_time>-0.1)
  }
  ddf <- ddf %>% mutate(p_level_fdr = as.factor(case_when(
    padj_fdr_term > .05 ~ '1',
    padj_fdr_term < .05 & padj_fdr_term > .01 ~ '2',
    padj_fdr_term < .01 & padj_fdr_term > .001 ~ '3',
    padj_fdr_term <.001 ~ '4')))
  ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
  ddf$`p, FDR-corrected` = ddf$p_level_fdr
  terms <- unique(ddf$term[ddf$effect=="fixed"])
  ddf$t <- ddf$evt_time
  ddf$t <- as.numeric(ddf$t)

  #define color schemes
  color_options <- c('#66CCEE','#c912ab','#e17014','#32af64','#644bc8','blue','firebrick1','darkseagreen1','gold')
  if (region_to_run == "subcortical") {
    color_scheme <- c(color_options[1:length(unique(ddf$network))])
    background_color='gray80'
    major_grid_color='gray90' 
    minor_grid_color='gray80' 
    hline_color='white' 
    vline_color='white'
  } else if (region_to_run == "cortical") {
    color_scheme <- c('#c912ab','#e17014','#644bc8')
    background_color='white' 
    major_grid_color='gray80'
    minor_grid_color='white'
    hline_color='gray40'
    vline_color='gray40'
  }
  
  #if (region_to_run == "subcortical") {if (split_hipp_amyg){color_scheme <- c('#66CCEE','#c912ab','#e17014','#32af64','#644bc8')
  #} else {color_scheme <- c('#32af64','#c912ab','#e17014','#644bc8')}
  #  background_color='gray80'
  #  major_grid_color='gray90' 
  #  minor_grid_color='gray80' 
  #  hline_color='white' 
  #  vline_color='white'
  #} else if (region_to_run == "cortical") {color_scheme <- c('#c912ab','#e17014','#644bc8') 
  #background_color='white' 
  #major_grid_color='gray80'
  #minor_grid_color='white'
  #hline_color='gray40'
  #vline_color='gray40'}
  
  #If want to flip statistic for PLOTS
  if (flip_stat) {ddf$statistic <- -ddf$statistic
  ddf$estimate <- -ddf$estimate}
  
  ########## Create the plots ##########
  for (fe in terms) {
    message(paste0("\nPlotting ",fe,"\n"))
    edf <- ddf %>% filter(term == paste(fe) & ddf$effect=='fixed')
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
    if (region_to_run == "subcortical") {
      if (flip_stat) {fname = paste0(fe,'-Plot-',structure,'-',model_name,"-FLIPPED.pdf")}
      else {fname = paste0(align_to,'-',fe,'-Plot-',structure,'-',model_name,".pdf")}
    } else if (region_to_run == "cortical") {
      if (flip_stat) {fname = paste0('vmPFC-',align_to,'-',fe,'-Plot-',structure,'-',model_name,"-FLIPPED.pdf")}
      else {fname = paste0('vmPFC-',align_to,'-',fe,'-Plot-',structure,'-',model_name,".pdf")}}
    if (exclude_prefeedback) {
      if (region_to_run == "cortical") {pdf(fname, width = 4.5, height = 4)
      } else {pdf(fname, width = 5, height = 4)}
    } else {pdf(fname, width = 6, height = 4)}
    pd <- position_dodge(0.2)
    gg<-ggplot(edf, aes(x=t, y=.data[[toplot]], color=network, group=network)) + #changed to flip y for groups compared to attempters 5/16 x=t, y=statistic, color=network)) + 
      geom_vline(xintercept = 0, lty = 'solid', color = vline_color, size = 1.5)+ xlab(epoch_label) + ylab('') +
      geom_hline(yintercept = 0, lty = 'solid',color = hline_color, size=1.5) + 
      geom_point(aes(size=`p, FDR-corrected`, alpha=`p, FDR-corrected`), position=pd) + #show.legend = FALSE
      scale_colour_manual(values = color_scheme) + 
      geom_line(aes(group=network, alpha = `p, FDR-corrected`), size=1.5, position=pd, show.legend = FALSE) + 
      geom_linerange(aes(ymin=estimate-std.error, ymax=estimate+std.error, alpha = `p, FDR-corrected`),show.legend = FALSE, position=pd) +
      scale_alpha_discrete(range=c(0.4, 1), drop=FALSE) + scale_size_discrete(range=c(2,6), drop=FALSE) # I think adding this last bit fixed my issue of the points getting resized when all factors weren't in the data
    gg <- gg  + theme_bw(base_size=13) +
      xlab(paste("time with respect to",align_to," (seconds)")) +
      ylab(ylabel) + 
      theme(legend.title = element_blank(),
            panel.grid.major = element_line(colour = major_grid_color), 
            panel.grid.minor = element_line(colour = minor_grid_color), 
            panel.background = element_rect(fill = background_color),
            axis.title.y = element_text(margin=margin(r=6)),
            axis.title.x = element_text(margin=margin(t=6)))
    print(gg)
    dev.off()
  }
}















