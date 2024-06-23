rm(list = ls())
library(ggplot2)
library(lme4)
library(lmerTest)
library(dplyr)
library(reshape)
library(nlme)
library(permuco)
library(ggseg)
library(stringr)
library(mediation)

code_location <- dirname(sys.frame(1)$ofile)
source(paste0(code_location, '/model_funs.R'))
source(paste0(code_location, '/graph_labels_colors.R'))
source(paste0(code_location, '/gen_heatmap.R'))
source(paste0(code_location, '/compare_models_fun.R'))
source(paste0(code_location, '/MSU_functions.R')) # from https://github.com/knickmeyer-lab/ORIGINs_ICV-and-Subcortical-volume-development-in-early-childhood
source(paste0(code_location, '/long_model_report_stats.R'))
source(paste0(code_location, '/diff_node_clust.R'))
source(paste0(code_location, '/brain_region_ggseg.R'))
source(paste0(code_location, '/long_graph.R'))
source(paste0(code_location, '/long_descriptive_stats.R'))
source(paste0(code_location, '/long_mediate.R'))



dir_in <- '<insert_input_directory>' # location of 'reorg' directory, containing reorganized MRI data and non-MRI data (i.e., 'other_dat.csv')
models <- c('lin', 'log') # model functions: lin, log, quad, asym
measures <- c('volume', 'wm') # volume, wm, area, thickness, meancurv, dti_fa, dti_md
trks <- c('ARC_L', 'SLF_L', 'ILF_L') # tracts if analyzing diffusion
covariates <- c('Sex', 'MEd', 'HLE', 'FHD') # covariates 
rating_t1 <- 1.5 # qc ratings for structure 
rating_dwi <- 0.5 # qc ratings for diffusion

              

# set additional parameters
trk_sz <- 'average' # degree to which tract is segmented: nodes, quarters, average'
hem <- 'lh' # lh, rh
dV <- 0.1 # inter-timepoint change threshold for preschool (PRE) to 2nd-3rd grade (REA) timepoints
ctr_str <- 'none' # control options: none, euler, eTIV 
use_inf <- 1 # require every participant to have at least one infant (INF) or toddler (TOD) observation: 0, 1
use3 <- 0 # require >2 time points: 0, 1, obsolete
vqc_enigma <- 0 # remove ids and timepoints based on visual qc with ENIGMA protocol: 0, 1
rating_vqc <- 1 # set rating value for vqc_enigma


regions <- c('bankssts', 'fusiform', 'inferiorparietal', 
             'middletemporal', 'parsopercularis', 
             'parstriangularis', 'superiortemporal', 'supramarginal') 

beh_names <- c('pre_wjiv_phonproc_standard',	'pre_wjiv_orallang_standard')
med_names <- c('beg_wrmt3_wrd_id_stnd',	'beg_towre2_composite_stnd')



# load all other data
other_dat <- read.csv( paste0( dir_in, '/reorg/other_dat.csv' ))



# run models
cat("\014")
for (meas in measures){
  cat('\n\n\n', paste0('---- Fitting growth curves and performing brain-behavior relationships on...', meas, '----'), '\n\n')
  i <- 0
  j <- 1
  while (i == 0){
    i = i + 1
    # rating based on measures
    if (meas == 'dti_fa' | meas == 'dti_md'){
      rating <- rating_dwi
      if (j < length(trks)){ # needed to loop through tracks
        i <- i - 1 # reset i to maintain while loop
      }
    } else {
      rating <- rating_t1
    }
    
    trk <- trks[j] # changes if dti_fa or dti_md is specified. ineffectual for other measures 
    j <- j + 1
    
    for (model in models){
      if (model == 'log'){
        log_model_fun(meas, hem, trk, trk_sz, dir_in, other_dat, covariates, regions, 
                         beh_names, rating, dV, ctr_str, use_inf, use3, vqc_enigma, rating_vqc)
      } else if (model == 'quad'){
        quad_model_fun(meas, hem, trk, trk_sz, dir_in, other_dat, covariates, regions,
                          beh_names, rating, dV, ctr_str, use_inf, use3, vqc_enigma, rating_vqc)
      } else if (model == 'lin'){
        lin_model_fun(meas, hem, trk, trk_sz, dir_in, other_dat, covariates, regions,
                         beh_names, rating, dV, ctr_str, use_inf, use3, vqc_enigma, rating_vqc)
      } else if (model == 'asym'){
        asym_model_fun(meas, hem, trk, trk_sz, dir_in, other_dat, covariates, regions,
                          beh_names, rating, dV, ctr_str, use_inf, use3, vqc_enigma, rating_vqc)
      } else {
        print('no model selected')
      }
   
      
      # generate heatmaps for covariates for lme models
      if (model == 'log' | model == 'quad' | model == 'lin'){
        gen_heatmap(stats_sum, meas, dir_in, regions, covariates)
      }
    }
    
    if (length(models) > 1){
      compare_models_fun(meas, hem, trk, trk_sz, dir_in, regions)
      next
    }
    
    # compute stats with multiple comparison corrections
    if (model != 'asym'){
    long_model_report_stats(meas, hem, trk_sz, dir_in, regions,
                                beh_names, covariates, stats_sum)
      if ((meas == 'dti_fa' | meas == 'dti_md') & trk_sz == 'nodes') {
        diff_node_clust(diff_clus_elem, stats_sum, meas, hem, trk, dir_in, beh_names)
      }
    } else {
    long_model_report_stats_asym(meas, hem, trk_sz, dir_in, regions,
                                beh_names, covariates, stats_sum)
      if ((meas == 'dti_fa' | meas == 'dti_md') & trk_sz == 'nodes') {
        diff_node_clust_asym(diff_clus_elem, meas, hem, trk, dir_in, beh_names)
      }
      next
    }
    
    # generate ggseg pdfs
    if (meas != 'dti_fa' & meas != 'dti_md'){
      if (length(indx_fdr) != 0){
        brain_region_ggseg(meas, hem, dir_in, indx_fdr)
      }
    }

    # run mediation
    if (model != 'asym'){
      long_mediate(med_elem, stats_sum, indx_fdr_bin, meas, hem, trk, trk_sz, dir_in, 
                   other_dat, regions, beh_names, med_names)
    }
  }
}


# -- Run after all desired measures have been run above --
# generate longitudinal plot
long_graph(dir_in, measures, hem, trks, covariates)

# generate descriptive statistics
long_descriptive_stats(dir_in, measures, hem, trks, covariates, beh_names)


