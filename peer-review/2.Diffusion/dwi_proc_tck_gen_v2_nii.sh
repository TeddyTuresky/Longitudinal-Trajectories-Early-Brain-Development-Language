#!/bin/bash

# Early Development Diffusion Pipeline
# This is the nifti-input analogue of dwi_proc_tck_gen_v2.sh

# trap keyboard interrupt (control-c)
trap control_c SIGINT

function Help {
    cat <<HELP
Usage:
`basename $0` -m MainDiffusionSequence -r ReversePhaseEncodingSequence -t NumberSimultSliceAcq -o OutputDirectory -a InfantFreeSurferDirectory -x BvecsFile -l BvalsFile

Compulsory arguments:
     -m:  Main sequence nifti file
     -r:  Directory containing raw dicoms with phase-encoding direction opposite main sequence
     -t:  Specify number of simultaneous slice acquisitions (1|2)
     -a:  Directory containing FreeSurfer style structure
     -o:  OutputDirectory: Full path to output directory
     -x:  Bvecs file for main sequence
     -l:  Bvals file for main sequence
Optional arguments:
     -f:  Method to generate FODs (msmt|ss3t)
     -p:  Abort after preprocessing 
--------------------------------------------------------------------------------------
script by T. Turesky, 2022
--------------------------------------------------------------------------------------

HELP
    exit 1
}


# Provide output for Help
if [[ "$1" == "-h" || $# -eq 0 ]]; then
    Help >&2
fi


# Read command line arguments
while getopts "h:m:r:t:o:a:x:l:f:p" OPT
  do
  case $OPT in
      h) #help
   Help
   exit 0
   ;;
      m)  # main sequence nifti
   main=$OPTARG
   ;;
      r)  # reverse pe directory
   rev=$OPTARG
   ;;
      t)  # specify number of simultaneous slices acquired
   tim=$OPTARG
   ;;
      a)  # anatomy directory
   anat=$OPTARG
   ;;
      o)  # output directory
   out=$OPTARG
   ;;
      x)  # bvecs file
   bvecs=$OPTARG
   ;;     
      l)  # bvals file
   bvals=$OPTARG
   ;;
      f)  # use (standard) msmt or ss3t
   fod=$OPTARG
   ;;
      p)  # only preprocess
   pre_only=1
   ;;
     \?) # getopts issues an error message
   echo "$USAGE" >&2
   exit 1
   ;;
  esac
done


# Ensure inputs exist
if [[ ! -f "${main}" ]]; then
    echo "Error: main sequence '${main}' does not exist."
    exit 1
fi

if [[ ! -f "${bvecs}" ]]; then
    echo "Error: bvecs file for main sequence '${bvecs}' does not exist."
    exit 1
fi

if [[ ! -f "${bvals}" ]]; then
    echo "Error: bvals file for main sequence '${bvals}' does not exist."
    exit 1
fi

if [[ ! -d "${rev}" ]]; then
    echo "Warning: reverse PE sequence '${rev}' does not exist or not selected. Not using topup."
    top=0
else
    top=1
fi

if [[ -d "${out}" ]]; then
     if [ "$(ls -A $out)" ]; then
	echo "Warning: chosen output directory not empty."
     fi
    # echo "Output directory '${out}' does not exist. Creating..." # below
fi

if [ "${tim}" != 1 ] && [ "${tim}" != 2 ]; then
    echo "Error: number of slices acquired simultaneously not specified or a number other than 1 or 2. Cannot run eddy with slice-to-volume correction."
    exit 1
fi

if [[ ! -d "${anat}" ]]; then
    echo "Error: directory containing anatomical images '${anat}' does not exist."
    exit 1
elif [[ ! -d "${anat}"/mri ]]; then
    echo "Error: directory '${anat}' exists, but does not contain FreeSurfer structure."
    exit 1
fi


# Set environmental variables - differs by choice of FOD tool
if [[ -v fod ]]; then

    case $fod in
        msmt)
        echo "You have opted for FOD generation with msmt."
        fd=1 # indicator for function below
        ;;

        ss3t)
        echo "You have opted for FOD generation with ss3t."
        fd=2 # indicator for function below
        ;;

        *)
        echo "Invalid specification for FOD generation. Defaulting to msmt."
        fd=1 # indicator for function below        
	;;
    esac

else
    echo "FOD generation tool not specified. Defaulting to msmt."
    fd=1 # indicator for function below
fi


# Generate output directory structure
pdir=${out}/preproc
mkdir -p $pdir


# Pull anatomical data and FreeSurfer outputs
mrconvert ${anat}/mri/orig.mgz ${pdir}/t1.nii.gz -datatype uint16 -force # maybe change to rawavg.mgz?
mrconvert ${anat}/mri/aseg.presurf.mgz ${pdir}/aseg.nii.gz -datatype uint16 -force # using aseg or aseg.presurf?
mrconvert ${anat}/mri/brainmask.mgz ${pdir}/brain.nii.gz -datatype uint16 -force 

fslorient -copyqform2sform ${pdir}/t1.nii.gz
fslorient -copyqform2sform ${pdir}/brain.nii.gz
fslorient -copyqform2sform ${pdir}/aseg.nii.gz

# Begin to preprocess diffusion data
mrconvert ${main} ${pdir}/dwi.mif -fslgrad ${bvecs} ${bvals} -strides -1,-2,+3,+4 -force 
#dwigradcheck ${main} -fslgrad ${bvecs} ${bvals} -export_grad_mrtrix ${pdir}/grad.txt
#mrconvert ${main} -grad ${pdir]/grad.txt ${pdir}/dwi.mif
dwidenoise ${pdir}/dwi.mif ${pdir}/den_dwi.mif -force


# Identify FSL-style slice-timing file to use based on number of simultaneous slice acquisitions and number of slices and continue preprocessing depending also on whether there is a reverse phase encoded sequence
n_slc=`mrinfo -size ${pdir}/dwi.mif | awk '{ print $3 }'`
slspec=`dirname $(which dwi_proc_tck_gen_v2.sh)`/slspecs/slspec_${tim}_${n_slc}.txt
echo Using ${slspec} as slice order

if [[ ! -f "${slspec}" ]]; then
    echo "Error: FSL slice order table '${slspec}', a prerequisite for eddy with slice-to-volume correction, does not exist."
    exit 1
else
    if [ ${top} = 1 ]; then
	echo performing susceptibility distortion correction
	mrconvert ${rev} ${pdir}/PA.mif -strides -1,-2,+3 -force

	# ensuring same number of slices in main and rev sequences
	n_slc_rev=`mrinfo -size ${pdir}/PA.mif | awk '{ print $3 }'`
	if [ "${n_slc}" != "${n_slc_rev}" ]; then
	    echo main and rev sequences have different numbers of slices. regridding rev image... 
	    mrgrid ${pdir}/PA.mif regrid -template ${pdir}/dwi.mif ${pdir}/PA2.mif -interp nearest -force # cannot do a direct overwrite so need the following two lines
	    rm ${pdir}/PA.mif 
	    cp ${pdir}/PA2.mif ${pdir}/PA.mif
	fi

	mrconvert ${pdir}/den_dwi.mif ${pdir}/1_b0.mif -coord 3 0 -axes 0,1,2 -force
	mrcat ${pdir}/1_b0.mif ${pdir}/PA.mif ${pdir}/1_b0_PA.mif -axis 3 -force
	dwifslpreproc ${pdir}/den_dwi.mif ${pdir}/pre_den_dwi.mif -pe_dir AP -rpe_pair -se_epi ${pdir}/1_b0_PA.mif -align_seepi -eddy_options " --mporder=6 --slm=linear --ol_nstd=5" -eddy_slspec=${slspec} -force
    else
	echo performing eddy without susceptibility distortion correction
	dwifslpreproc ${pdir}/den_dwi.mif ${pdir}/pre_den_dwi.mif -pe_dir AP -rpe_none -eddy_options " --mporder=6 --slm=linear --ol_nstd=5" -eddy_slspec=${slspec} -force
    fi
fi


# Finish preprocessing diffusion data
dwibiascorrect ants ${pdir}/pre_den_dwi.mif ${pdir}/N4_pre_den_dwi.mif -force


# Align T1 to DWI data
dwiextract ${pdir}/N4_pre_den_dwi.mif - -bzero | mrmath - mean ${pdir}/mean_bzero.mif -axis 3 -force
mrconvert ${pdir}/mean_bzero.mif ${pdir}/mean_bzero.nii.gz -force
antsRegistrationSyN.sh -d 3 -f ${pdir}/mean_bzero.nii.gz -m ${pdir}/brain.nii.gz -t r -o ${pdir}/brain2dwi
antsApplyTransforms -d 3 -i ${pdir}/brain.nii.gz -r ${pdir}/brain.nii.gz -o ${pdir}/brain2dwi.nii.gz -t ${pdir}/brain2dwi0GenericAffine.mat
antsApplyTransforms -d 3 -i ${pdir}/t1.nii.gz -r ${pdir}/t1.nii.gz -o ${pdir}/t12dwi.nii.gz -t ${pdir}/brain2dwi0GenericAffine.mat
antsApplyTransforms -d 3 -i ${pdir}/aseg.nii.gz -r ${pdir}/brain.nii.gz -o ${pdir}/aseg2dwi.nii.gz -t ${pdir}/brain2dwi0GenericAffine.mat -n NearestNeighbor
cp ${pdir}/N4_pre_den_dwi.mif ${out}/reg_N4_pre_den_dwi.mif


# Export nifti format with bvec and bval files for pyAFQ and generate a mask file
mrconvert ${out}/reg_N4_pre_den_dwi.mif ${out}/reg_N4_pre_den_dwi.nii.gz -export_grad_fsl ${out}/reg_N4_pre_den_dwi.bvec ${out}/reg_N4_pre_den_dwi.bval -strides -1,+2,+3,+4 -force
fslmaths ${pdir}/brain2dwi.nii.gz -bin ${pdir}/mask.nii.gz

# Abort after preprocessing if desired
echo ${pre_only}
if [ ${pre_only} = 1 ]; then
    exit 1
fi


# Generate 5tt and gmwmi files for ACT
5ttgen freesurfer ${pdir}/aseg2dwi.nii.gz ${out}/5tt.mif -nocrop -force
5tt2gmwmi ${out}/5tt.mif ${out}/5tt_gmwmi.mif -force


# Generate FODs and run tractography
echo "Generating FODs and running tractography."
pre=${out}/reg_N4_pre_den_dwi
dwi2mask ${pre}.mif ${pre}_mask.mif -force
dwi2response dhollander ${pre}.mif ${pre}_wm_resp.txt ${pre}_gm_resp.txt ${pre}_csf_resp.txt -fa 0.1 -lmax 0,8 -force

if [ ${fd} = 1 ]; then
    dwi2fod msmt_csd ${pre}.mif ${pre}_wm_resp.txt ${pre}_wm_fod.mif ${pre}_csf_resp.txt ${pre}_csf_fod.mif -mask ${pre}_mask.mif -force
elif [ ${fd} = 2 ]; then
    ${SS3T}/ss3t_csd_beta1 ${pre}.mif ${pre}_wm_resp.txt ${pre}_wm_fod.mif ${pre}_csf_resp.txt ${pre}_csf_fod.mif -mask ${pre}_mask.mif -force
else
    echo fod generation tool not specified
fi

mtnormalise ${pre}_wm_fod.mif ${pre}_wm_fod_norm.mif ${pre}_csf_fod.mif ${pre}_csf_fod_norm.mif -mask  ${pre}_mask.mif -force
tckgen ${pre}_wm_fod_norm.mif -algorithm iFOD1 -info -backtrack -crop_at_gmwmi -seed_gmwmi ${out}/5tt_gmwmi.mif -act ${out}/5tt.mif -select 3000000 ${out}/WholeBrainFG.tck -force
