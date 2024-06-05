#!/bin/python


import os
import os.path as op

import matplotlib.pyplot as plt
import nibabel as nib

from AFQ.api.group import GroupAFQ
import bids

bids_path='INSERT_SUB'


my_afq = GroupAFQ(
    bids_path,
    preproc_pipeline='mrtrix',
    import_tract={
        "suffix": "tractography",
        "scope": "my_tractography"
    })

my_afq.export_all()
