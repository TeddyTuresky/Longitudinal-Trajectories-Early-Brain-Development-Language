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
library(lm.beta)

code_location <- dirname(sys.frame(1)$ofile)
source(paste0(code_location, '/long_mod_gen.R'))
source(paste0(code_location, '/long_lmer_mods.R'))
source(paste0(code_location, '/asym_model_fun.R'))
source(paste0(code_location, '/mod_ages.R'))
source(paste0(code_location, '/mod_sages.R'))
source(paste0(code_location, '/graph_labels_colors.R'))
source(paste0(code_location, '/gen_heatmap.R'))
source(paste0(code_location, '/compare_models_bic_fun.R'))
source(paste0(code_location, '/compare_models_fun_covs.R'))
source(paste0(code_location, '/MSU_functions.R')) # from https://github.com/knickmeyer-lab/ORIGINs_ICV-and-Subcortical-volume-development-in-early-childhood
source(paste0(code_location, '/long_model_report_stats.R'))
source(paste0(code_location, '/diff_node_clust.R'))
source(paste0(code_location, '/brain_region_ggseg.R'))
source(paste0(code_location, '/long_mediate_covs.R'))
source(paste0(code_location, '/long_graph_wcal.R'))
source(paste0(code_location, '/graph_all_meas_regs.R'))
source(paste0(code_location, '/long_descriptive_stats.R'))
source(paste0(code_location, '/long_cov_beh_corr.R'))
source(paste0(code_location, '/gen_histograms_spec.R'))



dir_in <- '/Users/tht622/Library/CloudStorage/OneDrive-HarvardUniversity/hgse/dev_struct_2023.12.04' # '<location of "reorg" directory>' # should contain reorganized MRI data and non-MRI data (i.e., 'other_dat.csv')
models <- 'log' # model functions: lin, log, quad, asym
measures <- 'volume' # #volume, wm, area, thickness, meancurv, dti_fa, dti_md
trks <- c('ARC_L', 'SLF_L', 'ILF_L') # tracts if analyzing diffusion
covariates <- c('Sex', 'MEd', 'cohort') # covariates 'cohort', HLE', 'FHD'
covariates_interaction <- NULL # specify covariate(s) for interactions, otherwise set to NULL
covariates_mediation <- c('MEd', 'HLE', 'FHD') # covariates only for mediation
cov_var <- 'FHD' # covariate to test without, specify NA if not testing covariate significance
rating_t1 <- 1.5 # qc ratings for structure 
rating_dwi <- 0.5 # qc ratings for diffusion

              

# set additional parameters
trk_sz <- 'nodes' # degree to which tract is segmented: nodes, quarters, average'
hem <- 'lh' # lh, rh
dV <- 0.1 # annual change threshold for preschool (PRE) to 2nd-3rd grade (REA) timepoints
ctr_str <- 'none' # control options: none, euler, eTIV 
use_inf <- 0 # require every participant to have at least one infant (INF) or toddler (TOD) observation: 0, 1
use3 <- 0 # require >2 time points: 0, 1, obsolete
vqc_enigma <- 0 # remove ids and timepoints based on visual qc with ENIGMA protocol: 0, 1
rating_vqc <- 1 # set rating value for vqc_enigma


regions <- c('bankssts', 'fusiform', 'inferiorparietal', 
             'middletemporal', 'parsopercularis', 
             'parstriangularis', 'superiortemporal', 'supramarginal')



beh_names <- c('pre_phonproc_standard')
med_names <- c('beg_wrmt3_wrd_id_stnd',	'beg_wrmt3_wrdatt_stnd')


# load all other data
other_dat <- read.csv( paste0( dir_in, '/reorg/other_dat_tp.csv' ))



# run models
cat("\014")
for (meas in measures){
  cat('\n\n\n', paste0('---- Fitting growth curves and performing brain-behavior relationships on...', meas, ' ----'), '\n\n')
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
      if (model == 'lin' | model == 'log' | model == 'quad'){
        long_mod_gen(model, meas, hem, trk, trk_sz, dir_in, other_dat, covariates, 
                     covariates_interaction, cov_var, regions, beh_names, rating, 
                     dV, ctr_str, use_inf, use3, vqc_enigma, rating_vqc)
      } else if (model == 'asym'){
        asym_model_fun(meas, hem, trk, trk_sz, dir_in, other_dat, covariates, regions,
            beh_names, rating, dV, ctr_str, use_inf, use3, vqc_enigma, rating_vqc)
      }
   
      
      # generate heatmaps for covariates for lme models
      if (model != 'asym'){
        gen_heatmap(stats_sum, meas, dir_in, regions, covariates)
      }
      
      if (any(grepl(cov_var, covariates)) == 1){
        compare_models_fun_covs(model, meas, hem, trk, trk_sz, dir_in, cov_var, regions)
      }
    }
    
    if (length(models) > 1){
      compare_models_bic_fun(meas, hem, trk, trk_sz, dir_in, regions)
      next
    }
    
    
    
    # compute stats with multiple comparison corrections
    if (model != 'asym'){
    long_model_report_stats(meas, hem, trk_sz, dir_in, regions,
                                beh_names, covariates, covariates_interaction, stats_sum)
      if ((meas == 'dti_fa' | meas == 'dti_md') & trk_sz == 'nodes') {
        diff_node_clust(diff_clus_elem, stats_sum, meas, hem, trk, dir_in, beh_names)
      }
    } else {
    long_model_report_stats_asym(meas, hem, trk_sz, dir_in, regions,
                                beh_names, covariates, stats_sum)
      if ((meas == 'dti_fa' | meas == 'dti_md') & trk_sz == 'nodes') {
        diff_node_clust_asym(diff_clus_elem, meas, hem, trk, dir_in, beh_names)
      }
      #next
    }
    
    # run mediation
    if (model != 'asym'){
      long_mediate(med_elem, stats_sum, indx_fdr_bin, meas, hem, trk, trk_sz, dir_in, 
                   other_dat, covariates_mediation, regions, beh_names, med_names, ctr_str)
    }
    
    # generate histograms for random effects of brain trajectories
    if (model != 'asym'){
      gen_histograms_spec_brain_refs(dir_in, meas, hem, trk, med_elem, regions, beh_names)
    }
    
    # generate ggseg pdfs
    if (model != 'asym' & meas != 'dti_fa' & meas != 'dti_md'){
      
      # generate heatmap by cortical region
      brain_region_ggseg(meas, hem, dir_in, indx_est, 0)
      
      if (length(indx_fdr) != 0){
        brain_region_ggseg(meas, hem, dir_in, indx_fdr, 1)
      }
      if (model != 'asym'){
        if (length(indx_med_all) != 0){
          brain_region_ggseg_med(meas, hem, dir_in, indx_med_all)
        }
      }
    }

  }
}


# -- Run after all desired measures have been run above --
if (model != 'asym' & trk_sz == 'average' ){
  # generate all measures change graph (will use whatever model was last run, except asymptotic), need to select average trk_sz
  graph_all_meas_regs(dir_in, measures, hem, trk)
  
  # generate longitudinal plot
  long_graph_wcal(dir_in, measures, hem, trks, covariates, use_inf)
  
  # generate descriptive statistics
  long_descriptive_stats(dir_in, measures, hem, trks, covariates, beh_names, med_names, use_inf)
  
  # test for correlations among behavioral measures using final dataset
  long_cov_beh_corr(dir_in, hem, covariates_mediation, beh_names, med_names)
  
  # generate histogram for key behavioral variables
  gen_histograms_spec_beh(dir_in, hem, other_dat)
  
}