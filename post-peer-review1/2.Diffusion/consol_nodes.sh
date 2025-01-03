#!/bin/bash

# consolidates diffusion estimates in preparation for statistical analysis

# Input:
#  1 - full path of py(Baby)AFQ directory



dir_in = ${1}

e=(`echo ${dir_in}/*`); 
mkdir -p ${dir_in}/nodes

for i in ${e[@]}; do 
	echo $i
	
	cp ${i}/derivatives/afq/sub-01/ses-01/sub-01_ses-01_dwi_space-RASMM_model-probCSD_algo-AFQ_desc-profiles_dwi.csv ${dir_in}/nodes/${i}_nodes.csv 
done