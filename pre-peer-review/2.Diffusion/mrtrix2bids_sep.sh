#!/bin/bash

# Generates BIDS structure for input to pyAFQ from mrtrix/freesurfer preprocessed DWI and tractography files

# Example: mrtrix2bids_sep.sh </path2/mrtrix-freesurfer_folder> </path2/pyAFQ_folder>

fun=`dirname $(which mrtrix2bids_sep.sh)`

cp -R ${fun}/temp ${2}

cp ${1}/reg_N4_pre_den_dwi.nii.gz ${2}/derivatives/mrtrix/sub-01/ses-01/dwi/sub-01_ses-01_dwi.nii.gz
cp ${1}/reg_N4_pre_den_dwi.bvec ${2}/derivatives/mrtrix/sub-01/ses-01/dwi/sub-01_ses-01_dwi.bvec
cp ${1}/reg_N4_pre_den_dwi.bval ${2}/derivatives/mrtrix/sub-01/ses-01/dwi/sub-01_ses-01_dwi.bval

cp ${1}/WholeBrainFG.tck ${2}/derivatives/my_tractography/sub-01/ses-01/dwi/sub-01_ses-01-dwi_tractography.tck
