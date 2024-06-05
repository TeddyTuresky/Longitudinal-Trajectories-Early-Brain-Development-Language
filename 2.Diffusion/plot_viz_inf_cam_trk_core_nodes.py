# the code below is adapted from https://yeatmanlab.github.io/pyAFQ/ 

import os
import os.path as op
import nibabel as nib
import numpy as np
import sys
import csv

from dipy.io.streamline import load_trk
from dipy.tracking.streamline import transform_streamlines

from fury import actor, window
from fury.colormap import create_colormap


dir_in = sys.argv[1]
sub_id = sys.argv[2]
dir_out = sys.argv[3]
tract = sys.argv[4]
hem = sys.argv[5]
node_file = sys.argv[6]

node_pre = os.path.splitext(node_file)[0]
nodes = np.genfromtxt(node_file, delimiter=',')


color_trk = np.array((0, 0.8, 0.4))


study_path = os.path.join(dir_in, sub_id)

deriv_path = op.join(
    study_path, "derivatives")

afq_path = op.join(
    deriv_path,
    'afq',
    'sub-01',
    'ses-01')

bundle_path = op.join(afq_path,
                      'clean_bundles')


fa_img = nib.load(op.join(afq_path,
'sub-01_ses-01_dwi_model-DTI_desc-FA_dwi.nii.gz'))
fa = fa_img.get_fdata()
sft_trk = load_trk(op.join(bundle_path,
     'sub-01_ses-01_dwi_space-RASMM_model-probCSD_algo-AFQ_desc-'+tract+hem+'_tractography.trk'), fa_img)



t1w_img = nib.load(op.join(afq_path,'sub-01_ses-01_dwi_desc-b0_dwi.nii.gz'))
t1w = t1w_img.get_fdata()
sft_trk.to_rasmm()
trk_t1w = transform_streamlines(sft_trk.streamlines,
                                np.linalg.inv(t1w_img.affine))


def lines_as_tubes(sl, line_width, **kwargs):
    line_actor = actor.line(sl, **kwargs)
    line_actor.GetProperty().SetRenderLinesAsTubes(1)
    line_actor.GetProperty().SetLineWidth(line_width)
    return line_actor


color_bck = (0.99, 0.99, 0.99) # (0.75, 0.75, 0.75) 
trk_actor = lines_as_tubes(trk_t1w, 2, colors=color_trk, opacity=0.1)



from dipy.tracking.streamline import set_number_of_points
core_trk = np.median(np.asarray(set_number_of_points(trk_t1w, 100)), axis=0)


from dipy.stats.analysis import afq_profile
sft_trk.to_vox()
trk_profile = afq_profile(fa, sft_trk.streamlines, affine=np.eye(4))



core_color_bck = np.tile(np.array(color_bck), (100,1))


for i in range(len(nodes)):
	core_color_bck[int(nodes[i]+4)] = color_trk # offset because we removed first (and last) 5 nodes and -1 converting numbering to python



core_trk_actor = lines_as_tubes(
    [core_trk],
    200,
    colors = core_color_bck 
    )



def slice_volume(data, x=None, y=None, z=None):
    slicer_actors = []
    slicer_actor_z = actor.slicer(data)
    if x is not None:
        slicer_actor_x = slicer_actor_z.copy()
        slicer_actor_x.display_extent(
            x, x,
            0, data.shape[1] - 1,
            0, data.shape[2] - 1)
        slicer_actors.append(slicer_actor_x)
    return slicer_actors


slicers = slice_volume(t1w, x=t1w.shape[0]//2, z=t1w.shape[-1]//3)


scene = window.Scene()

scene.add(core_trk_actor)


if hem == 'L':
    scene.set_camera(position=(250, 30, 40),
                 focal_point=(55, 45, 30),
                 view_up=(0, 0, 1))
elif hem == 'R':
    scene.set_camera(position=(-140, 43, 40),
                 focal_point=(55, 45, 30),
                 view_up=(0, 0, 1))
else:
    print('error: hemisphere not set appropriately')

#scene.zoom(0.25)

scene.background((255, 255, 255))

window.record(scene, out_path=os.path.join(dir_out, sub_id+'_'+tract+hem+node_pre+'.png'), size=(2400, 2400))
