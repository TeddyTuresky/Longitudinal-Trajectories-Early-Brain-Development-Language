#!/bin/bash

# Finishes FS recon-all -autorecon1 and -autorecon2-wm pipelines with some adjustments from FS 7.3


fp=${1}
ifp=${2}
sub=${3}




cd ${fp}/mri

# Resumes -autorecon1 pipeline using iFS files (incomplete because we already have many iFS files that can be used instead)
mri_nu_correct.mni --i orig.mgz --o nu.mgz --uchar transforms/talairach.xfm --n 2 --ants-n4
mri_normalize -g 1 -seed 1234 -mprage nu.mgz T1.mgz
mri_mask ${fp}/mri/T1.mgz ${ifp}/mri/brainmask.mgz ${fp}/mri/brainmask.mgz
#mri_normalize -g 1 -seed 1234 -mprage nu.mgz T1.mgz
#mri_em_register -skull nu.mgz /Users/tht622/freesurfer/average/RB_all_withskull_2020_01_02.gca transforms/talairach_with_skull.lta
#mri_watershed -T1 -brain_atlas /Users/tht622/freesurfer/average/RB_all_withskull_2020_01_02.gca transforms/talairach_with_skull.lta T1.mgz brainmask.auto.mgz
mri_em_register -uns 3 -mask brainmask.mgz nu.mgz ${FREESURFER_HOME}/average/RB_all_2020-01-02.gca transforms/talairach.lta
mri_ca_normalize -c ctrl_pts.mgz -mask brainmask.mgz nu.mgz ${FREESURFER_HOME}/average/RB_all_2020-01-02.gca transforms/talairach.lta norm.mgz
mri_normalize -seed 1234 -mprage -aseg aseg.presurf.mgz -mask brainmask.mgz norm.mgz brain.mgz


# Performs most of the functions in -autorecon2
mri_mask -T 5 brain.mgz brainmask.mgz brain.finalsurfs.mgz 
mri_fill -a ../scripts/ponscc.cut.log -xform transforms/talairach.lta -segmentation aseg.presurf.mgz wm.mgz filled.mgz  -ctab ${FREESURFER_HOME}/SubCorticalMassLUT.txt 
cp filled.mgz filled.auto.mgz
mri_pretess filled.mgz 255 norm.mgz filled-pretess255.mgz
mri_tessellate filled-pretess255.mgz 255 ../surf/lh.orig.nofix
rm -f filled-pretess255.mgz
mri_pretess filled.mgz 127 norm.mgz filled-pretess127.mgz
mri_tessellate filled-pretess127.mgz 127 ../surf/rh.orig.nofix
rm -f filled-pretess127.mgz

cd ${fp}/surf

mris_extract_main_component lh.orig.nofix lh.orig.nofix
mris_extract_main_component rh.orig.nofix rh.orig.nofix
mris_smooth -nw -seed 1234 lh.orig.nofix lh.smoothwm.nofix
mris_smooth -nw -seed 1234 rh.orig.nofix rh.smoothwm.nofix
mris_inflate -no-save-sulc lh.smoothwm.nofix lh.inflated.nofix
mris_inflate -no-save-sulc rh.smoothwm.nofix rh.inflated.nofix
mris_sphere -q -p 6 -a 128 -seed 1234 -in 3000 lh.inflated.nofix lh.qsphere.nofix # added -in 3000 from iFS
mris_sphere -q -p 6 -a 128 -seed 1234 -in 3000 rh.inflated.nofix rh.qsphere.nofix # added -in 3000 from iFS

cp lh.orig.nofix lh.orig
cp rh.orig.nofix rh.orig
cp lh.qsphere.nofix lh.qsphere
cp rh.qsphere.nofix rh.qsphere

mris_topo_fixer -mgz -warnings ${sub} lh
mris_topo_fixer -mgz -warnings ${sub} rh

mv lh.orig_corrected lh.orig.premesh
mv rh.orig_corrected rh.orig.premesh
rm lh.orig rh.orig lh.orig_corrected rh.orig_corrected

mris_euler_number lh.orig.premesh
mris_euler_number rh.orig.premesh
#defect2seg --s ${sub} --cortex
mris_remesh --remesh --iters 3 --input lh.orig.premesh --output lh.orig
mris_remesh --remesh --iters 3 --input rh.orig.premesh --output rh.orig
mris_remove_intersection lh.orig lh.orig
mris_remove_intersection rh.orig rh.orig
mris_autodet_gwstats --o autodet.gw.stats.lh.dat --i ../mri/brain.finalsurfs.mgz --wm ../mri/wm.mgz --surf lh.orig.premesh
mris_autodet_gwstats --o autodet.gw.stats.rh.dat --i ../mri/brain.finalsurfs.mgz --wm ../mri/wm.mgz --surf rh.orig.premesh

mris_make_surfaces -output .preaparc -soap -orig_white orig -aseg aseg.presurf -cover_seg ${fp}/mri/aseg.presurf.mgz -noaparc -whiteonly -mgz -T1 brain.finalsurfs ${sub} lh
mris_make_surfaces -output .preaparc -soap -orig_white orig -aseg aseg.presurf -cover_seg ${fp}/mri/aseg.presurf.mgz -noaparc -whiteonly -mgz -T1 brain.finalsurfs ${sub} rh

mri_label2label --label-cortex lh.white.preaparc ../mri/aseg.presurf.mgz 0 ../label/lh.cortex.label
mri_label2label --label-cortex lh.white.preaparc ../mri/aseg.presurf.mgz 1 ../label/lh.cortex+hipamyg.label
mri_label2label --label-cortex rh.white.preaparc ../mri/aseg.presurf.mgz 0 ../label/rh.cortex.label
mri_label2label --label-cortex rh.white.preaparc ../mri/aseg.presurf.mgz 1 ../label/rh.cortex+hipamyg.label
mris_smooth -n 3 -nw -seed 1234 lh.white.preaparc lh.smoothwm
mris_smooth -n 3 -nw -seed 1234 rh.white.preaparc rh.smoothwm
mris_inflate lh.smoothwm lh.inflated
mris_inflate rh.smoothwm rh.inflated
mris_curvature -w -seed 1234 lh.white.preaparc
ln -s lh.white.preaparc.H lh.white.H
ln -s lh.white.preaparc.K lh.white.K
mris_curvature -seed 1234 -thresh .999 -n -a 5 -w -distances 10 10 lh.inflated
mris_curvature -w -seed 1234 rh.white.preaparc
ln -s rh.white.preaparc.H rh.white.H
ln -s rh.white.preaparc.K rh.white.K
mris_curvature -seed 1234 -thresh .999 -n -a 5 -w -distances 10 10 rh.inflated
