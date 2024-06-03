# earlyLongitudinalBrainDevelopmentLanguage

This repository houses code (or links to code) used for the following study:

    Turesky et al. Longitudinal trajectories of brain development from infancy to school-age and their relationship to language skill

Broadly, the study involved structural and diffusion processing pipelines followed by statistical analyses. An inventories of the code used for each pipeline are outlined below:

Structure:


    .
    ├── diffusion_pipeline_gen.sh                     <-- generates tract reconstructions from raw diffusion images using MRTrix, VistaSoft, and AFQ
        ├── vista_preprocessing.m                     <-- aligns diffusion and T1 images
        ├── mrtrix2babyAFQ.m                          <-- runs AFQ fiber segmentation
            ├── dti_end_tract.m                       <-- converts MRtrix .tck file to WholeBrainFG.mat for AFQ fiber segmentation
    ├── babyAFQ_est_vis.m                          <-- generates tract profiles with diffusion estimates and visualizations for figures
    ├── nonpar_boot_sp_corr.m                         <-- runs semipartial correlations (code fragments taken from https://github.com/yeatmanlab/AFQ) 
    ├── AFQ_MultiCompCorrectionSemiPartSpearman.m     <-- runs multiple comparison corrections for brain-behavior relations, adjusted from https://github.com/yeatmanlab/AFQ
    ├── hleMedsOpen.R                                 <-- tests for mediation



Diffusion:


The scripts enclosed in this repository rely on the following packages:
. Mrtrix3 | Mrtrix3Tissue

. Matlab R2021a and the following packages (the following order of priority helps):
. . Vistasoft


Dependences for optional calls
. acpcdetect


Adjustments to LD_LIBRARY_PATH
. GCC libraries


Statistical analyses:
    .
    ├── bb_struct_reorg.R                             <-- reformats structural estimates from FreeSurfer 
    ├── bb_dwi_reorg_r4.R                             <-- reformats diffusion estimates from py(Baby)AFQ
    ├── longitudinal_model_control.R                  <-- performs longitudinal modelling, including data cleaning, curve fitting and visualizations, and associations with outcomes
        ├── model_funs.R                              <-- runs longitudinal models: linear, logarithmic, quadratic, or asymptotic
        ├── graph_labels_colors.R                     <-- specifies colors used in graphs
        ├── gen_heatmap.R                             <-- runs heatmaps depicting covariate contributions to models
        ├── compare_models_fun.R                      <-- compares fits among longitudinal models (does not include asymptotic functions)
        ├── MSU_longitudinal_models_functions.R       <-- leverages a subset of functions provided here: https://github.com/knickmeyer-lab/ORIGINs_ICV-and-Subcortical-volume-development-in-early-childhood (asymptotic models only)
        ├── long_model_report_stats_fun               <-- generates and reports statistics corrected for multiple comparisons (for structural and average- or quarter-based diffusion analyses)  
        ├── diff_node_clust_r5.R                      <-- generates and reports statistics corrected for multiple comparisons (for node-based diffusion analyses)
        ├── brain_region_ggseg_fun_r3.R               <-- depicts significant brain regions (for structure only)
        ├── long_graph_r2.R                           <-- generates graph showing longitudinal sample
        ├── long_descriptive_stats_r6.R               <-- reports descriptive stats for study


required packages: dplyr, reshape, stringr, ggplot2, ggseg, lme4, lmerTest, nlme, permuco





