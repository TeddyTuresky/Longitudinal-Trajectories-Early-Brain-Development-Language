#!/bin/bash

# Finishes FS recon-all -autorecon3 with some adjustments from FS 7.3


fp=${1}
fun=${2}
sub=${3}



recon-all -autorecon3 -subjid ${sub} -expert ${fun}/expert.opts # use expert options to cause recon-all to crash at Pial Surface step
mris_make_surfaces -grad_dir 1 -intensity .3 -output .tmp -pial_offset .25 -nowhite -noaparc -aseg aseg.presurf -cover_seg ${fp}/mri/aseg.presurf.mgz -orig_pial white ${sub} lh # use call from iFS
mris_make_surfaces -grad_dir 1 -intensity .3 -output .tmp -pial_offset .25 -nowhite -noaparc -aseg aseg.presurf -cover_seg ${fp}/mri/aseg.presurf.mgz -orig_pial white ${sub} rh # use call from iFS
mv ${fp}/surf/lh.pial.tmp ${fp}/surf/lh.pial.T1
mv ${fp}/surf/rh.pial.tmp ${fp}/surf/rh.pial.T1
ln -sf ${fp}/surf/lh.pial.T1 ${fp}/surf/lh.pial
ln -sf ${fp}/surf/rh.pial.T1 ${fp}/surf/rh.pial
recon-all -autorecon3-T2pial -noT2pial -subjid ${sub} # run the rest of the recon-all -autorecon3 pipeline from the T2pial step but without the T2pial step - limited options for intermediate pipelines
