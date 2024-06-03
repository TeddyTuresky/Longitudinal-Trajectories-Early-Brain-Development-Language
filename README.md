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







