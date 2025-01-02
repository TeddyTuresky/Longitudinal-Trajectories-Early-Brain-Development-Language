#!/bin/bash
module load cuda/10.2.89-fasrc01 freesurfer/6.0.0-fasrc01 ANTs/2.3.1-fasrc01 Anaconda3/5.0.1-fasrc01 netpbm/10.73.20-fasrc01 ncf dcm2niix/2019_09_04-ncf fsl/6.0.4-ncf
# Anaconda3/2022.05 matlab/R2021a-fasrc01 

export tools=/n/gaab_mri_l3/Lab/DMC-Gaab2/tools/tkt_tools/babyBOLD_dmri_proc_pyafq
export FSLDIR=/n/home_fasse/tturesky/usr/local/fsl
export GCC=/n/sw/eb/apps/centos7/GCCcore/11.2.0/lib64
#export ARTHOME=${tools}/art
export SS3T=${tools}/mrtrix_3tissue/MRtrix3Tissue/bin
export PATH=${tools}:/n/home_fasse/tturesky/.conda/pkgs/mrtrix3-3.0.3-h2bc3f7f_0/bin:$FSLDIR/bin:$PATH  # $ARTHOME/bin:
export LD_LIBRARY_PATH=$GCC:/n/helmod/apps/centos7/Core/Anaconda3/5.0.1-fasrc01/x/lib:$LD_LIBRARY_PATH # /n/sw/eb/apps/centos7/Anaconda3/2022.05/lib64

source ${FSLDIR}/etc/fslconf/fsl.sh


echo Done Loading Dependencies
