import numpy as np
import pandas as pd
import nibabel as nib
import os

#Input
dseg_T1w_subfields = snakemake.input[0]

#Params
voxel_shape = snakemake.params[0]
#Wildcards
subj = snakemake.wildcards.subject
hemi = snakemake.wildcards.hemi

#Output
output_volumetry_T1w_subfields_tsv = str(snakemake.output[0])

#Separate subfields into separate variables, count number of voxels, convert to volumes mm3
subfields = ['Subiculum', 'CA1', 'CA2', 'CA3', 'CA4', 'DG']
subfield_nvox = {}
subfield_volume = {}
total_hipp_volume = 0

#Load in dseg_subfields
dseg_subfields = nib.load(dseg_T1w_subfields).get_fdata()

for sf_idx, sf_name in enumerate(subfields):
    sf_nvox = np.count_nonzero(dseg_subfields == sf_idx+1)
    subfield_nvox[sf_name] = sf_nvox
    subfield_volume[sf_name] = subfield_nvox[sf_name] * voxel_shape[0] * voxel_shape[1] * voxel_shape[2]
    total_hipp_volume += subfield_volume[sf_name]
subfield_volume['total_hipp_volume'] = total_hipp_volume


volumes_df = pd.DataFrame.from_dict(subfield_volume, orient='index').T
volumes_df = volumes_df.rename(columns={'Subiculum': 'Subiculum_volume', 'CA1': 'CA1_volume', 'CA2': 'CA2_volume', 'CA3': 'CA3_volume', 'CA4': 'CA4_volume', 'DG': 'DG_volume', 'total_hipp_volume': 'total_hipp_volume'})


#Write to tsv
volumes_df.to_csv(output_volumetry_T1w_subfields_tsv, sep='\t', index=False)
