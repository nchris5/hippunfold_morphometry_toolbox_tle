#Group Volumetry

import numpy as np
import pandas as pd
import re
import os
import csv

#Inputs
L_volumetry_T1w_subfields_tsv = snakemake.input.L_volumetry_T1w_subfields_tsv
R_volumetry_T1w_subfields_tsv = snakemake.input.R_volumetry_T1w_subfields_tsv

#Params
demographics_tsv = pd.read_csv(snakemake.params[0], sep='\t')
participant_ids = demographics_tsv['participant_id'].values.tolist()
tle_lateralization = demographics_tsv['tle_lateralization'].values.tolist() #Control R L BL
mri_normal = demographics_tsv['MRI_normal'].values.tolist() #N ABN

#Output
group_volumetry_T1w_subfields_tsv_dir = str(snakemake.output[0])


#Define function to loop through the expanded input from snakemake to create a grouped dataframe of all subjects distance for L and R separately
#  This loops through input files and adds a row for each subject, skipping the header row
def subj_to_group_df(tsvs, hemi):
    group_hemi = []
    for i, file in enumerate(tsvs):
        subj_df = pd.read_table(file, sep='\t')
        print(subj_df)
        group_hemi.append({'participant_id': participant_ids[i], 'tle_lateralization': tle_lateralization[i], 'MRI_normal': mri_normal[i], 'hemi': hemi, 'Subiculum_volume': subj_df['Subiculum_volume'].values[0], 'CA1_volume': subj_df['CA1_volume'].values[0], 'CA2_volume': subj_df['CA2_volume'].values[0], 'CA3_volume': subj_df['CA3_volume'].values[0], 'CA4_volume': subj_df['CA4_volume'].values[0], 'DG_volume': subj_df['DG_volume'].values[0], 'total_hipp_volume': subj_df['total_hipp_volume'].values[0]})
    group_hemi_df = pd.DataFrame(group_hemi)
    return group_hemi_df

group_volumetry_L = subj_to_group_df(L_volumetry_T1w_subfields_tsv, 'L')
group_volumetry_R = subj_to_group_df(R_volumetry_T1w_subfields_tsv, 'R')
group_volumetry_LR = pd.concat([group_volumetry_L, group_volumetry_R]).reset_index(drop=True)


#Define function to separate grouped input into ipsilateral_TLE, contralateral_TLE, bilateral_TLE, and Controls
#For each of these groups, Group by participant_id to extract all the volumetry values for each unique participant_id, so that a new df is created where each row contains a list of all the values for the volumetry along with the unique participant_id
def subgroup_split_volumetry(group_v_df_LR):
    subgroups_dict_dfs = {}
    ipsilateral = group_v_df_LR.query('tle_lateralization in ("R", "L") & hemi == tle_lateralization').copy()
    contralateral = group_v_df_LR.query('tle_lateralization in ("R", "L") & hemi != tle_lateralization').copy()
    bilateral = group_v_df_LR.query('tle_lateralization == "BL"').copy()
    allAffected = pd.concat([ipsilateral, bilateral]).reset_index(drop=True)
    all_TLE = pd.concat([ipsilateral, contralateral, bilateral]).reset_index(drop=True)
    controls = group_v_df_LR.query('tle_lateralization == "Control"').copy()
    subgroups_dict_dfs['ipsilateral_TLE'] = ipsilateral
    subgroups_dict_dfs['contralateral_TLE'] = contralateral
    subgroups_dict_dfs['bilateral_TLE'] = bilateral
    subgroups_dict_dfs['allAffected_TLE'] = allAffected
    subgroups_dict_dfs['all_TLE'] = all_TLE
    subgroups_dict_dfs['controls'] = controls
    return subgroups_dict_dfs

subgroups_volumetry_LR_dict_dfs = subgroup_split_volumetry(group_volumetry_LR)

#Write each of the subgroups' volumetry dfs to separate tsvs
if not os.path.exists(group_volumetry_T1w_subfields_tsv_dir):
    os.mkdir(group_volumetry_T1w_subfields_tsv_dir)

for key, df in subgroups_volumetry_LR_dict_dfs.items():
    df.to_csv(group_volumetry_T1w_subfields_tsv_dir+'/group_'+key+'_volumetry.tsv', sep='\t', index=False)


