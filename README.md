# earlyLongitudinalBrainDevelopmentLanguage

This repository houses code (or links to code) used for the following study:

    Turesky et al. Longitudinal trajectories of brain development from infancy to school-age and their relationship to language skill

Broadly, the study involved structural and diffusion processing pipelines followed by statistical analyses. An inventory of the code used for each pipeline is provided below:

### 1. Structure:

    .
    ├── reFS.sh                                 <-- runs infant brain morphometry pipeline
        ├── iFS_wrap.sh                         <-- sets environmental variables for iFS and runs it *this file will need adjustments based on the user's compute setup
        ├── ibeat2aseg.m                        <-- merges iFS and iBEATv2.0 tissue segmentations
            ├── aseg_labels2coords.m            <-- converts aseg labels to coordinates
        ├── aseg2wm.m                           <-- generates FS white matter file using iBEATv2.0 segmentation
        ├── fs_autorecon2_end.sh                <-- runs FS autorecon2 with adjustments
        ├── fs_autorecon3_wrap.sh               <-- runs FS autorecon3 with adjustments, uses expert.opts file included
    ├── consol_stats.sh                         <-- consolidates structural estimates from FS 
    ├── struct_reorg.R                          <-- reorganizes structural estimates

Dependencies: iBEATv2.0 Docker, Infant FreeSurfer, FreeSurfer 7.3, Matlab  
Additional requirements:The folder containing the above scripts also need to be added to your shell PATH, iBEATv2.0 (https://github.com/iBEAT-V2/iBEAT-V2.0-Docker) segmentations need to be generated beforehand.  


### 2. Diffusion:

    .
    ├── dwi_proc_tck_gen_v2.sh                     <-- preprocesses data and runs fiber tracking, *uses eddy_cuda, which requires an slspec file be created in advance, please see FSL documentation
    ├── mrtrix2bids_sep.sh                         <-- reorganizes directory structure to prepare for py(Baby)AFQ
    ├── baby-pyafq-gen.py                          <-- runs tract segmentation for pyBabyAFQ, bids_path needs to be entered in file
    ├── pyafq-gen.py                               <-- runs tract segmenation for pyAFQ, bids_path needs to be entered in file
    ├── plot_viz_inf_cam_trk.py                    <-- visualizes tracts; code adapted from https://yeatmanlab.github.io/pyAFQ/
    ├── plot_viz_inf_cam_trk_core_nodes.py         <-- visualizes tract cores with significant nodes, code adapted from https://yeatmanlab.github.io/pyAFQ
    ├── consol_nodes.sh                            <-- consolidates diffusion estimates from py(Baby)AFQ
    ├── dwi_reorg.R                                <-- reorganizes diffusion estimates

Dependences: MRtrix3, FSL, ANTs, FreeSurfer, GCC libraries for LD_LIBRARY_PATH, pyAFQ  
Additional requirements: the folder containing the above scripts needs to be added to your shell PATH, our first script (indirectly) calls the cuda implementation of eddy, which also requires access to a GPU (we used Nvidia A100), FreeSurfer-style T1 segmentation need to be generated beforehand.


### 3. Statistics:

    .
    ├── longitudinal_model_control.R            <-- set paramaters for statistical analyses visualizations
        ├── model_funs.R                        <-- runs longitudinal models - linear, logarithmic, quadratic, or asymptotic - including data cleaning, curve fitting and visualizations, and associations with outcomes
        ├── graph_labels_colors.R               <-- specifies colors used in graphs
        ├── gen_heatmap.R                       <-- runs heatmaps depicting covariate contributions to models
        ├── compare_models_fun.R                <-- compares fits among longitudinal models, does not include asymptotic functions
        ├── MSU_functions.R                     <-- leverages a subset of functions provided here: https://github.com/knickmeyer-lab/ORIGINs_ICV-and-Subcortical-volume-development-in-early-childhood (asymptotic models only)
        ├── long_model_report_stats.R           <-- generates and reports statistics corrected for multiple comparisons, for structural and average- or quarter-based diffusion analyses
        ├── diff_node_clust.R                   <-- generates and reports statistics corrected for multiple comparisons, for node-based diffusion analyses
        ├── brain_region_ggseg.R                <-- depicts significant brain regions, for structure only
        ├── long_graph.R                        <-- generates graph showing longitudinal sample
        ├── long_descriptive_stats.R            <-- reports descriptive stats for study


Required packages: dplyr, reshape, stringr, ggplot2, ggseg, lme4, lmerTest, nlme, permuco  
Additional requirements: asymptotic functions also require code supplied here: https://github.com/knickmeyer-lab/ORIGINs_ICV-and-Subcortical-volume-development-in-early-childhood.

