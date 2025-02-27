#!/bin/bash

# This script runs an infant brain morphometry pipeline using iBEATv2 and standard FreeSurfer (FS). The following pipeline assumes that iBEATv2 tissue segmentations have already been generated.
#   1. Run FS up to step 15 (https://surfer.nmr.mgh.harvard.edu/fswiki/recon-all#StepDescriptionSummaries)
#   2. Merge FS aseg.nii.gz with iBEATv2 T1-iso-skullstripped-tissue.nii.gz using the Matlab scripts:
#      a. ibeat2aseg.m - uses tissue labels from iBEATv2 to relabel FS aseg.mgz gray/white matter except for subcortical regions and then writes as aseg.presurf.mgz file in FS mri folder
#      b. aseg2wm.m - generates wm.mgz files for input to FS
#   3. Finish FS recon-all pipeline with autorecon2-wm and autorecon3



# Inputs
#   1. subjid
#   2. iBEATv2 tissue segmentation (1=csf, 2=gm, 3=wm)




#---------------------------------------------------------------------------------------------------------------------------


sub=$1
fp=${SUBJECTS_DIR}/${sub}

m_fp=\'${fp}/mri\' # output location for ibeat2aseg.m
m_ib=\'${2}\'
m_aseg=\'${fp}/mri/aseg.presurf.nii\' # output from ibeat2aseg.m and input to aseg2wm.m
fun=`dirname $(which ibeat2aseg.m)`
m_fun=\'${fun}\'

echo Using ... $FREESURFER_HOME

# Check inputs
if [[ ! -f "${2}" ]]; then
    echo "iBEATv2 tissue segmentation not specified."
    exit 1
fi

if [ $# -ne 2 ]; then
    echo "Two arguments not specified. Participant ID or iBEATv2 filepath may be missing."
    exit 1
fi

# Step 1. 
echo Running first steps of recon-all
recon-all -all -subjid ${sub} -nonormalization2


# Step 2.
echo Merging iBEATv2 tissue segmentation with FS aseg.presurf and generating new aseg.presurf and wm files for FS
matlab -nodesktop -nosplash -r "addpath(${m_fun}); ibeat2aseg(${m_ib}, ${m_fp}); aseg2wm(${m_aseg}); exit;"


# Step 3.
echo Finishing FS recon-all pipeline
mri_normalize -seed 1234 -mprage -aseg ${fp}/mri/aseg.presurf.mgz -mask ${fp}/mri/brainmask.mgz ${fp}/mri/norm.mgz ${fp}/mri/brain.mgz
recon-all -autorecon2-wm -autorecon3 -subjid ${sub}
