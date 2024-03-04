#AI - plot_mixed_by_simple_AI.R; uses ddf generated from Build_mixed_by_AI.R
#Last updated Feb 7, 2024

rm(list=ls())

library(wesanderson)
library(pracma)
library(ggplot2)
library(tidyverse)
library(viridis)
library(grid)
library(viridis) 

#Choose flags:
toalign <- "clock" #options are clock or response
toprocess <- "network" # Options I want are atlas_value (every ROI and hemisphere separately), subregion (collapse across subregions and hemispheres), network (collapse across networks and hemispheres), structure (coarser groups, combined across hemispheres), structure_hipp (with anteiror and posterior hippocampus separated)
colorby <- "structure" #options are subregion and structure
toplot <- "estimate" #options: are estimate and statistic
model <- "model_5" # The name of the model you want to plot (name defined in "Build_mixed_by.R" script)
region_to_run = "subcortical" #options: cortical (vmPFC), subcortical 
flip_stat = 0 #to flip statistic for plots (i.e. if model set up to look at attempters vs. depressed but you want to flip to see depressed vs. attempters)
split_hipp_amyg = TRUE #if want to separate out hippocampus and amygdala
exclude_prefeedback = 0 #set to 1 if you want to exclude the pre-feedback results from plots

#Set the working directory based on the options
subcortical_cache_dir = '/Users/angela/Documents/RESEARCH/ANALYSES/Medusa_analysis/subcortical'  # base directory, where MEDUSA output files are
cortical_cache_dir = '/Users/angela/Documents/RESEARCH/ANALYSES/Medusa_analysis/cortical'
vmpfc_hc_cache_dir = '/Users/angela/Documents/RESEARCH/ANALYSES/Medusa_analysis/vmPFC_hipp_interactions'

#Define where data is located and where plots will be saved
if(region_to_run == "subcortical") {
  if (toprocess == "structure") {dir <- paste0(subcortical_cache_dir,'/Structure_CombinedHemis')
  dir.create(file.path(subcortical_cache_dir, '/Structure_CombinedHemis'), showWarnings = FALSE)
  } else if (toprocess == "subregion") {dir <- paste0(subcortical_cache_dir,'/Subregion_CombinedHemis')
  dir.create(file.path(subcortical_cache_dir, '/Subregion_CombinedHemis'), showWarnings = FALSE)
  } else if (toprocess == "network") {dir <- paste0(subcortical_cache_dir,'/Network_CombinedHemis')
  dir.create(file.path(subcortical_cache_dir, '/Network_CombinedHemis'), showWarnings = FALSE)
  } else if (toprocess == "atlas_value") {dir <- paste0(subcortical_cache_dir,'/Network_CombinedHemis')
  dir.create(file.path(subcortical_cache_dir, '/Network_CombinedHemis'), showWarnings = FALSE)}
} else if(region_to_run == "cortical") {
  if (toprocess == "structure") {dir <- paste0(cortical_cache_dir,'/Structure_CombinedHemis')
  dir.create(file.path(cortical_cache_dir, '/Structure_CombinedHemis'), showWarnings = FALSE)
  } else if (toprocess == "subregion") {dir <- paste0(cortical_cache_dir,'/Subregion_CombinedHemis')
  dir.create(file.path(cortical_cache_dir, '/Subregion_CombinedHemis'), showWarnings = FALSE)
  } else if (toprocess == "network") {dir <- paste0(cortical_cache_dir,'/Network_CombinedHemis')
  dir.create(file.path(cortical_cache_dir, '/Network_CombinedHemis'), showWarnings = FALSE)
  } else if (toprocess == "atlas_value") {dir <- paste0(cortical_cache_dir,'/Network_CombinedHemis')
  dir.create(file.path(cortical_cache_dir, '/Network_CombinedHemis'), showWarnings = FALSE)}
} else if(toprocess == "vmpfc_hc") {dir <- vmpfc_hc_cache_dir}
setwd(dir)

#Load the Rdata file with mixed_by results
load(dir(full.names=T,pattern=(paste("*-",toalign,"-",model,".Rdata", sep=""))))

if(toalign == "response"){epoch_label <- 'response'
} else if(toalign == "clock"){epoch_label <- 'clock'}

if(toplot == "statistic") {ylabel <- 't-statistic'
}else if(toplot == "estimate") {ylabel <- 'Coefficient (AU)'}

#If you're plotting the individual ROIs, load and join results with the labels file to make it easier to interpret
#Note: ROIs 202 and 204 are both labeled as R BLA in the labels file so plots are a little wonky because of that
labels <- read_csv("~/Documents/RESEARCH/ANALYSES/fMRI/Dec2022/extracted_values/12Dec2022/region_labels_244.csv") %>% mutate(roi_num = as.numeric(roi_num)) %>% inner_join(read_csv("~/Documents/RESEARCH/ANALYSES/fMRI/Dec2022/extracted_values/12Dec2022/region_lookup_244.csv"), by = "roi_num")
labels <- labels %>% mutate(atlas_value = as.numeric(roi_num))
if (toprocess == "atlas_value") {ddf$coef_df_reml <- ddf$coef_df_reml %>% inner_join(labels, by="atlas_value")
} else if (toprocess == "subregion") {ddf$coef_df_reml <- ddf$coef_df_reml %>% inner_join(labels, by="subregion")
} else if (toprocess == "structure") {ddf$coef_df_reml <- ddf$coef_df_reml[,-7]
} else if (toprocess == "structure_hipp") {ddf$coef_df_reml <- ddf$coef_df_reml[,-7]
} 

# ddf is the output structure; reml and ml are two different convergence methods (use reml - contains AIC and BIC andn log likelihood)
# ddf$coeff_df_reml - has all fixed and random effects; value is the estimate and estimated error is std.error
#consolidate the ddf to just the coef_df_reml
ddf <- ddf$coef_df_reml

#ddf$atlas_value <- as.factor(ddf$atlas_value)
#ddf$subregion <- as.factor(ddf$subregion)
ddf$network <- as.factor(ddf$network)

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
if (region_to_run == "subcortical") {if (split_hipp_amyg){color_scheme <- c('#66CCEE','#c912ab','#e17014','#32af64','#644bc8')
} else {color_scheme <- c('#32af64','#c912ab','#e17014','#644bc8')}
  background_color='gray80'
  major_grid_color='gray90' 
  minor_grid_color='gray80' 
  hline_color='white' 
  vline_color='white'
} else if (region_to_run == "cortical") {color_scheme <- c('#c912ab','#e17014','#644bc8') 
background_color='white' 
major_grid_color='gray80'
minor_grid_color='white'
hline_color='gray40'
vline_color='gray40'}

#If want to flip statistic for PLOTS
if (flip_stat) {ddf$statistic <- -ddf$statistic
ddf$estimate <- -ddf$estimate}

for (fe in terms) {
  message(paste0("\nPlotting ",fe,"\n"))
  edf <- ddf %>% filter(term == paste(fe) & ddf$effect=='fixed')
  termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
  if (region_to_run == "subcortical") {if (flip_stat) {fname = paste0(toalign,'-',fe,'-Plot-',toprocess,'-',model,'-Colorby-',colorby,"FLIPPED.pdf")}
    else {fname = paste0(toalign,'-',fe,'-Plot-',toprocess,'-',model,'-Colorby-',colorby,".pdf")}
  } else if (region_to_run == "cortical") {if (flip_stat) {fname = paste0('vmPFC-',toalign,'-',fe,'-Plot-',toprocess,'-',model,'-Colorby-',colorby,"FLIPPED.pdf")}
    else {fname = paste0('vmPFC-',toalign,'-',fe,'-Plot-',toprocess,'-',model,'-Colorby-',colorby,".pdf")}}
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
    xlab(paste("time with respect to",toalign," (seconds)")) +
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
