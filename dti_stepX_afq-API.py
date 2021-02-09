# %%
# import os
# import os.path as op
# import matplotlib.pyplot as plt
# import nibabel as nib
from AFQ import api
import bids

# %%
bids_path = "/Users/nmuncy/Projects/emu_diff/data"
layout = bids.BIDSLayout(bids_path, derivatives=True)
print(layout)
print(layout.derivatives)

# %%
# check bids layout
my_vistasoft = layout.derivatives["vistasoft"]
vistasoft_files = my_vistasoft.get(extension=".nii.gz")
vistasoft_file = vistasoft_files[0]
print(vistasoft_file.get_entities())

# %%
# check bids truth
validator = bids.BIDSValidator()
vistasoft_full_path = vistasoft_file.path
vistasoft_relative_path = vistasoft_full_path.split(layout.root)[-1]
print(validator.is_bids(vistasoft_relative_path))

# %%
# prepare AFQ object
tract_list = ["SLF", "ARC", "CST", "FP"]
my_afq = api.AFQ(
    bids_path,
    dmriprep="vistasoft",
    bundle_info=tract_list,
    segmentation_params={"nb_streamlines": 1000},
)

# %%
my_afq.export_all()
