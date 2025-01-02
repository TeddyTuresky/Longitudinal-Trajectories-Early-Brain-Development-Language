#!/bin/bash

# Inputs:
#  1 - full path to subject list file
#  2 - output directory

# Check inputs
if [[ ! -f "${1}" ]]; then
    echo "Subject list file not specified or does not exist."
    exit 1
fi

if [[ ! -d "${2}" ]]; then
    echo "Output directory not specified. Saving to current directory"
    out_dir=./
else
    out_dir=${2}
fi




for hm in lh rh; do

  for m in volume area thickness thicknessstd meancurv gauscurv foldind curvind; do

    aparcstats2table --subjectsfile=${1} --hemi $hm --meas $m --tablefile ${out_dir}/${hm}_${m}_aparc_stats.txt

  done

done

asegstats2table --subjectsfile=${1} --meas volume --stats wmparc.stats --common-segs --tablefile ${out_dir}/volume_wm_stats.txt
