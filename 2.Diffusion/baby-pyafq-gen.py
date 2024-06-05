#!/bin/python


import os.path as op
import plotly

from AFQ.api.group import GroupAFQ
import AFQ.api.bundle_dict as abd
import AFQ.data.fetch as afd

from AFQ.definitions.image import RoiImage, ImageFile

bids_path='/n/holyscratch01/LABS/gaab_mri_l3/Lab/diff_long/baby-pyafq/INSERT_SUB'


my_afq = GroupAFQ(
    bids_path,
    preproc_pipeline='mrtrix',
    reg_template=afd.read_pediatric_templates(
    )["UNCNeo-withCerebellum-for-babyAFQ"],
    bundle_info=abd.PediatricBundleDict(),
    reg_subject_spec='b0',
    import_tract={
        "suffix": "tractography",
        "scope": "my_tractography"
    },
    segmentation_params={
        "filter_by_endpoints": False
    })

my_afq.export_all()


