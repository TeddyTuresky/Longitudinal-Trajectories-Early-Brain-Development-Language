#!/bin/bash

# This script runs an infant brain morphometry pipeline using iBEATv2, infant FreeSurfer (iFS), and standard FreeSurfer (FS). The following pipeline assumes that iBEATv2 tissue segmentations have already been generated.
#   1. Runs FS up to step 15 (https://surfer.nmr.mgh.harvard.edu/fswiki/recon-all#StepDescriptionSummaries)
#   2. Deletes all files in touch and script folders and wm.mgz and aseg.presurf.mgz files in mri folder
#   3. Converts FS orig.mgz images to mprage.nii.gz images and sets up for iFS 
#   4. Runs iFS in full
#   5. Merges iFS aseg.nii.gz with iBEATv2 T1-iso-skullstripped-tissue.nii.gz using the Matlab scripts:
#      a. ibeat2aseg.m - uses tissue labels from iBEATv2 to relabel iFS aseg.nii.gz gray/white matter except for subcortical regions and then writes as aseg.presurf.mgz file in FS mri folder
#      b. aseg2wm.m - generates wm.mgz files for input to FS
#   6. Resumes FS recon-all pipeline with some adjustments
#   7. Finishes FS recon-all pipeline with autorecon3



# Inputs
#   1. subjid
#   2. iBEATv2 tissue segmentation (1=csf, 2=gm, 3=wm)
#   3. age for infant FS




#---------------------------------------------------------------------------------------------------------------------------


sub=$1
fp=${SUBJECTS_DIR}/${sub}
if_dir=`dirname ${SUBJECTS_DIR}`/iFS
ifp=${if_dir}/${sub}

m_fp=\'${fp}/mri\' # output location for ibeat2aseg.m
m_ib=\'${2}\'
m_ifp=\'${ifp}\'
m_aseg=\'${fp}/mri/aseg.presurf.nii\' # output from ibeat2aseg.m and input to aseg2wm.m
fun=`dirname $(which ibeat2aseg.m)`
m_fun=\'${fun}\'

echo Using ... $FREESURFER_HOME

# Check inputs
if [[ ! -f "${2}" ]]; then
    echo "iBEATv2 tissue segmentation not specified."
    exit 1
fi

if [ $# -ne 3 ]; then
    echo "Three arguments not specified. Participant ID and/or age may be missing."
    exit 1
fi

# Step 1. 
echo Running first steps of recon-all
#recon-all -all -nofill -subjid ${sub} # the -nofill flag stops recon-all at step 15 (does not complete step 15)
recon-all -all -subjid ${sub} -nonuintensitycor


# Step 2.
echo Removing files based on initial FS run
# rm ${fp}/touch/* ${fp}/mri/wm*.mgz ${fp}/mri/aseg.presurf.mgz ${fp}/mri/brain.finalsurfs.mgz # ${fp}/scripts/* 
rm ${fp}/mri/transforms/* ${fp}/mri/orig_nu.mgz ${fp}/mri/mri_nu_correct.mni.log

# Step 3.
echo Setting up for infant FS
mkdir -p ${ifp}
mri_convert -i ${fp}/mri/orig.mgz -o ${ifp}/mprage.nii.gz


# Step 4.
echo Running infant FS
iFS_wrap.sh ${if_dir} $1 $3

echo Using ... $FREESURFER_HOME


# Step 5.
echo Merging iBEATv2 tissue segmentation with infant FS aseg and generating aseg.presurf and wm files for FS
# printf "${m_ib}\n${m_ifp}\n${m_fp}\n${m_aseg}\n"
matlab -nodesktop -nosplash -r "addpath(${m_fun}); ibeat2aseg(${m_ib}, ${m_ifp}, ${m_fp}); aseg2wm(${m_aseg}); exit;"


# Step 6.
echo Finishing FS recon-all -autorecon2-wm pipeline, including going back and performing some -autorecon1 steps using iFS files
cp ${ifp}/mri/transforms/talairach*xfm ${fp}/mri/transforms
fs_autorecon2_end.sh ${fp} ${ifp} ${sub}


# Step 7.
echo Running FS recon-all -autorecon3 pipeline with adjustments
fs_autorecon3_wrap.sh ${fp} ${fun} ${sub}
