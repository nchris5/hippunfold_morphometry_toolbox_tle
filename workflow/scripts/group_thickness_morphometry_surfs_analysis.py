#Group thickness

import numpy as np
import pandas as pd
import re
import os
import csv
import subprocess

#Inputs
###Thickness values
L_tsv_subfield_thickness_spaceUnfold_mPD_atEachAP = snakemake.params[1]
R_tsv_subfield_thickness_spaceUnfold_mPD_atEachAP = snakemake.params[4]

###FFT thickness values
L_tsv_FFT_subfield_thickness_spaceUnfold_mPD_atEachAP = snakemake.params[2]
R_tsv_FFT_subfield_thickness_spaceUnfold_mPD_atEachAP = snakemake.params[5]

###PSD of FFT thickness values
L_tsv_PSDofFFT_subfield_thickness_spaceUnfold_mPD_atEachAP = snakemake.params[3]
R_tsv_PSDofFFT_subfield_thickness_spaceUnfold_mPD_atEachAP = snakemake.params[6]

#Demographics table
demographics_tsv = pd.read_csv(snakemake.params[0], sep='\t')
participant_ids = demographics_tsv['participant_id'].values.tolist()
tle_lateralization = demographics_tsv['tle_lateralization'].values.tolist() #Control R L BL
mri_normal = demographics_tsv['MRI_normal'].values.tolist() #N ABN

#Wildcards: Surface
surface = 'midthickness' 

#Output dir
group_tsv_subfield_thickness_spaceUnfold_mPD_atEachAP_dir = str(snakemake.output[0])


#Define function to loop through the expanded input from snakemake to create a grouped dataframe of all subjects thickness for L and R separately
#  This loops through input files and adds a row for each subject
def subj_to_group_tsv_reader(tsvs, hemi, metric_name):
    group_hemi = []
    for i, file in enumerate(tsvs):
        with open(file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                metric_val = row[0].replace("'", "") # remove single quotes
                group_hemi.append({'participant_id': participant_ids[i], 'tle_lateralization': tle_lateralization[i], 'MRI_normal': mri_normal[i], 'hemi': hemi, metric_name: metric_val})
    group_hemi_df = pd.DataFrame(group_hemi)
    return group_hemi_df

#Define function to separate grouped input into ipsilateral_TLE, contralateral_TLE, bilateral_TLE, and Controls
#For each of these groups, Group by participant_id to extract all the mPDthickness_mAPcentered values for each unique participant_id, so that a new df is created where each row contains a list of all the values for the metric in a list along  with the unique participant_id
def subgroup_split(group_LR_df, metric_name):
    subgroups_dict_dfs = {}
    ipsilateral = group_LR_df.query('tle_lateralization in ("R", "L") & hemi == tle_lateralization').copy()
    ipsilateral = ipsilateral.groupby('participant_id')[metric_name].agg(list).reset_index()
    subgroups_dict_dfs['ipsilateral_TLE'] = ipsilateral
    contralateral = group_LR_df.query('tle_lateralization in ("R", "L") & hemi != tle_lateralization').copy()
    contralateral = contralateral.groupby('participant_id')[metric_name].agg(list).reset_index()
    subgroups_dict_dfs['contralateral_TLE'] = contralateral
    bilateral = group_LR_df.query('tle_lateralization == "BL"').copy()
    bilateral = bilateral.groupby(['participant_id', 'hemi'])[metric_name].agg(list).reset_index().drop('hemi', axis=1)
    subgroups_dict_dfs['bilateral_TLE'] = bilateral
    allAffected = pd.concat([ipsilateral, bilateral]).reset_index(drop=True)
    subgroups_dict_dfs['allAffected_TLE'] = allAffected
    all_TLE = pd.concat([ipsilateral, contralateral, bilateral]).reset_index(drop=True)
    subgroups_dict_dfs['all_TLE'] = all_TLE
    controls = group_LR_df.query('tle_lateralization == "Control"').copy()
    controls = controls.groupby(['participant_id', 'hemi'])[metric_name].agg(list).reset_index().drop('hemi', axis=1)
    subgroups_dict_dfs['controls'] = controls
    return subgroups_dict_dfs

###1) Group thickness values
group_mPDthickness_mAPcentered_L = subj_to_group_tsv_reader(L_tsv_subfield_thickness_spaceUnfold_mPD_atEachAP, 'L', 'mPDthickness_mAPcentered')
group_mPDthickness_mAPcentered_R = subj_to_group_tsv_reader(R_tsv_subfield_thickness_spaceUnfold_mPD_atEachAP, 'R', 'mPDthickness_mAPcentered')
#LR combined
group_mPDthickness_mAPcentered_LR = pd.concat([group_mPDthickness_mAPcentered_L, group_mPDthickness_mAPcentered_R]).reset_index(drop=True)
#Subgroups dictionary
subgroups_mPDthickness_mAPcentered_LR_dict_dfs = subgroup_split(group_mPDthickness_mAPcentered_LR, 'mPDthickness_mAPcentered')

###2) FFT of thickness values
group_FFT_mPDthickness_mAPcentered_L = subj_to_group_tsv_reader(L_tsv_FFT_subfield_thickness_spaceUnfold_mPD_atEachAP, 'L', 'FFT_mPDthickness_mAPcentered')
group_FFT_mPDthickness_mAPcentered_R = subj_to_group_tsv_reader(R_tsv_FFT_subfield_thickness_spaceUnfold_mPD_atEachAP, 'R', 'FFT_mPDthickness_mAPcentered')
#LR combined
group_FFT_mPDthickness_mAPcentered_LR = pd.concat([group_FFT_mPDthickness_mAPcentered_L, group_FFT_mPDthickness_mAPcentered_R]).reset_index(drop=True)
#Subgroups dictionary
subgroups_FFT_mPDthickness_mAPcentered_LR_dict_dfs = subgroup_split(group_FFT_mPDthickness_mAPcentered_LR, 'FFT_mPDthickness_mAPcentered')

###3) PSD of FFT of thickness values
group_PSDofFFT_mPDthickness_mAPcentered_L = subj_to_group_tsv_reader(L_tsv_PSDofFFT_subfield_thickness_spaceUnfold_mPD_atEachAP, 'L', 'PSDofFFT_mPDthickness_mAPcentered')
group_PSDofFFT_mPDthickness_mAPcentered_R = subj_to_group_tsv_reader(R_tsv_PSDofFFT_subfield_thickness_spaceUnfold_mPD_atEachAP, 'R', 'PSDofFFT_mPDthickness_mAPcentered')
#LR combined
group_PSDofFFT_mPDthickness_mAPcentered_LR = pd.concat([group_PSDofFFT_mPDthickness_mAPcentered_L, group_PSDofFFT_mPDthickness_mAPcentered_R]).reset_index(drop=True)
#Subgroups dictionary
subgroups_PSDofFFT_mPDthickness_mAPcentered_LR_dict_dfs = subgroup_split(group_PSDofFFT_mPDthickness_mAPcentered_LR, 'PSDofFFT_mPDthickness_mAPcentered')



#Define function to write each of the subgroups' dfs to separate tsvs
#Subgrouped tsvs have columns participant_id and metric_name (list of that metric's values with rows=#participants in that subgroup
#Bash subprocess is called here to remove single quotes around each element of the metric_name list (also removes whitespace between [ ( since this is an issue with the fft values when writing 
def save_subgroup_split_tsvs(outdir, subgroups_dict_dfs, metric_name):
    for key, df in subgroups_dict_dfs.items():
        filename = outdir+'/group_surf-'+surface+'_'+key+'_'+metric_name+'.tsv'
        df.to_csv(filename, sep='\t', index=False)
        subprocess.run(['sed', '-i', "-e", "s/'//g", "-e", "s/\[[[:space:]]/[/g", filename])

#Create output directory
if not os.path.exists(group_tsv_subfield_thickness_spaceUnfold_mPD_atEachAP_dir):
    os.mkdir(group_tsv_subfield_thickness_spaceUnfold_mPD_atEachAP_dir)

#Write subgroups to tsvs (with quotes removed and leading whitespace in lists removed)
save_subgroup_split_tsvs(group_tsv_subfield_thickness_spaceUnfold_mPD_atEachAP_dir, subgroups_mPDthickness_mAPcentered_LR_dict_dfs, 'mPDthickness_mAPcentered')
save_subgroup_split_tsvs(group_tsv_subfield_thickness_spaceUnfold_mPD_atEachAP_dir, subgroups_FFT_mPDthickness_mAPcentered_LR_dict_dfs, 'FFT_mPDthickness_mAPcentered')
save_subgroup_split_tsvs(group_tsv_subfield_thickness_spaceUnfold_mPD_atEachAP_dir, subgroups_PSDofFFT_mPDthickness_mAPcentered_LR_dict_dfs, 'PSDofFFT_mPDthickness_mAPcentered')

