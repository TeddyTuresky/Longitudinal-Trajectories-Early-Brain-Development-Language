#!/bin/bash


# change path parameters for iFS (and not FS)

module load ncf fsl/6.0.4-ncf
export FREESURFER_HOME=/n/gaab_mri_l3/Lab/DMC-Gaab2/tools/tkt_tools/infant_freesurfer/freesurfer
source $FREESURFER_HOME/SetUpFreeSurfer.sh
export FS_LICENSE=/n/gaab_mri_l3/Lab/DMC-Gaab2/tools/tkt_tools/infant_freesurfer/freesurfer/license.txt
export SUBJECTS_DIR=${1}


# Run infant FS
infant_recon_all --s ${2} --age ${3}
