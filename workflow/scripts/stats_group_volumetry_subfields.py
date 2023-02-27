#Volumetry statistical analyses

import numpy as np
import pandas as pd
from scipy.stats import f_oneway, ttest_ind #For ANOVA and t-test respectively
import matplotlib.pyplot as plt
import csv
import os

###Input dir
group_volumetry_T1w_subfields_tsv_dir = str(snakemake.input)

###Params 
#Subgroup names [ipsilateral_TLE, contralateral_TLE, bilateral_TLE, allAffected_TLE, all_TLE, controls]
group_names = snakemake.params[0]
#Subfield volumetry names [Subiculum_volume, CA1_volume, CA2_volume, CA3_volume, CA4_volume, DG_volume]
group_sf_volumetry_names = snakemake.params[1]

###Output dir
stats_group_volumetry_T1w_subfields_tsv_dir = str(snakemake.output)


#Define function to put all input subgroups volumetry data into a dictionary of dfs
def create_groups_dict_dfs_volumetry(groups_dir, names):
    groups_dict = {}
    for subgroup in names:
        f = groups_dir+'/group_'+subgroup+'_volumetry.tsv'
        groups_dict[subgroup] = pd.read_csv(f, sep='\t')
    return groups_dict
#create_groups_dict_dfs_volumetry (dictionary with keys=subgroup_names and values=df_subfields_volumes
groups_dict_dfs_volumetry = create_groups_dict_dfs_volumetry(group_volumetry_T1w_subfields_tsv_dir, group_names)

#Define function to calculate descriptive stats (mean, std, sterr, and var) for all subgroups and allsubfields; return these metrics in a nested dictionary subgroup[subfield[mean, std, sterr, var]]
#ddof is set to 1 as pandas defaults to 0 ddof

def calc_desc_stats_volumetry(groups_dict):
    groups_stats_dict = {}
    for subgroup, df in groups_dict.items():
        df_subfields = df.loc[:, ['Subiculum_volume', 'CA1_volume', 'CA2_volume', 'CA3_volume', 'CA4_volume', 'DG_volume']]
        subfields_stats_dict = {}
        for subfield in df_subfields.columns:
            subfield_mean = df_subfields[subfield].mean()
            subfield_std = df_subfields[subfield].std(ddof=1)
            subfield_stderr = subfield_std / np.sqrt(len(df))
            subfield_var = df_subfields[subfield].var(ddof=1)
            subfields_stats_dict[subfield] = {'mean': subfield_mean, 'std': subfield_std, 'stderr': subfield_stderr, 'variance': subfield_var}
        groups_stats_dict[subgroup] = subfields_stats_dict
    return groups_stats_dict
#calc_desc_stats_volumetry (nested dictionary with upperlevel_keys=subgroup_names, lowerlevel_keys=subfield_names, values=[mean, std, stderr, variance']
groups_subfields_nestedDict_dfs_stats_volumetry = calc_desc_stats_volumetry(groups_dict_dfs_volumetry)

#Function for downstream bonferroni correction during subfield volumetry statistical analyses
def bonferroni_correction(p_values, alpha):
    n_comparisons = len(p_values)
    corrected_alpha = alpha / n_comparisons
    corrected_p_values = [p_value * n_comparisons if p_value * n_comparisons < 1 else 1 for p_value in p_values]
    return corrected_p_values

#Define function for t-tests or ANOVA between subgroups whole hippocampus volumes
#  If subgroups(wholeHc) == 2 are provided --> t-test
#  If subgroups(wholeHc) > 2 are provided --> ANOVA
#  If subfields(Subiculum, CA1, CA2, CA3, CA4, DG) == True --> ONLY t-test followed by bonferroni correction
def ttest_ANOVA_volumetry_subgroups(*args, subfields=False):
    n_groups = len(args)
    if n_groups == 2 and not subfields: #t-test
        t, p = ttest_ind(*args)
        return t, p
    elif n_groups > 2 and not subfields: #ANOVA
        f, p = f_oneway(*args)
        return f, p
    elif n_groups == 2 and subfields: #Multiple subfields t-tests --> bonferroni correction
        t_p_subfields_dict = {}
        for subfield in group_sf_volumetry_names:
            groupA_sf = args[0][subfield]
            groupB_sf = args[1][subfield]
            t, p = ttest_ind(groupA_sf, groupB_sf)
            t_p_subfields_dict[subfield] = {'t': t, 'p': p}
        #Bonferroni correction returning p_corrected [1]
        p_subfields = [p_subfield['p'] for p_subfield in t_p_subfields_dict.values()]
        pCorr_subfields = bonferroni_correction(p_subfields, alpha=0.05)
        for i, subfield in enumerate(group_sf_volumetry_names):
            t_p_subfields_dict[subfield]['pCorr_bonferroni'] = pCorr_subfields[i]
        return t_p_subfields_dict

#Total Hc volume ANOVA: ipsilateral vs. contralateral vs. controls 
totalHc_volume_ipsVSconVScontrol_ANOVA_f, totalHc_volume_ipsVSconVScontrol_ANOVA_p = ttest_ANOVA_volumetry_subgroups(groups_dict_dfs_volumetry['ipsilateral_TLE']['total_hipp_volume'], groups_dict_dfs_volumetry['contralateral_TLE']['total_hipp_volume'], groups_dict_dfs_volumetry['controls']['total_hipp_volume'])
print('totalHc_volume_ipsVSconVScontrol_ANOVA_p = ', totalHc_volume_ipsVSconVScontrol_ANOVA_p,'\n')

#Subfields (Subiculum, CA1, CA2, CA3, CA4, DG) ttests and bonferroni corrected p-values: ipsilateral vs. contralateral
subfield_volumes_ipsVScon_ttestBonferroni_t_p_pCorr_dict = ttest_ANOVA_volumetry_subgroups(groups_dict_dfs_volumetry['ipsilateral_TLE'], groups_dict_dfs_volumetry['contralateral_TLE'], subfields=True)
print('subfield_volumes_ipsVScon_ttestBonferroni_pCorr_dict =', subfield_volumes_ipsVScon_ttestBonferroni_t_p_pCorr_dict, '\n')

##Subfields (Subiculum, CA1, CA2, CA3, CA4, DG) ttests and bonferroni corrected p-values: ipsilateral vs. controls
subfield_volumes_ipsVScontrols_ttestBonferroni_t_p_pCorr_dict = ttest_ANOVA_volumetry_subgroups(groups_dict_dfs_volumetry['ipsilateral_TLE'], groups_dict_dfs_volumetry['controls'], subfields=True)
print('subfield_volumes_ipsVScontrols_ttestBonferroni_t_p_pCorr_dict =', subfield_volumes_ipsVScontrols_ttestBonferroni_t_p_pCorr_dict, '\n')



#Create output dir and save stats tsvs
if not os.path.exists(stats_group_volumetry_T1w_subfields_tsv_dir):
    os.mkdir(stats_group_volumetry_T1w_subfields_tsv_dir)

#1) Save wholeHc ipsilateral vs. contralateral vs. controls ANOVA f and p
with open(stats_group_volumetry_T1w_subfields_tsv_dir+'/totalHc_volume_ipsilateral_vs_contralateral_vs_controls_ANOVA_f_and_p.tsv', 'w', newline='') as file:
    writer = csv.writer(file, delimiter='\t')
    writer.writerow([totalHc_volume_ipsVSconVScontrol_ANOVA_f, totalHc_volume_ipsVSconVScontrol_ANOVA_p])
#2) Save subfields ipsilateral vs. contralateral bonferroni-corrected t, p, pCorr
with open(stats_group_volumetry_T1w_subfields_tsv_dir+'/subfield_volumes_ipsilateral_vs_contralateral_ttestBonferroni_t_p_pCorr.tsv', 'w', newline='') as file:
    writer = csv.writer(file, delimiter='\t')
    for key, value in subfield_volumes_ipsVScon_ttestBonferroni_t_p_pCorr_dict.items():
        row = [key, value['t'], value['p'], value['pCorr_bonferroni']]
        writer.writerow(row)
#3) Save subfields ipsilateral vs. controls bonferroni-corrected t, p, pCorr
with open(stats_group_volumetry_T1w_subfields_tsv_dir+'/subfield_volumes_ipsVScontrols_ttestBonferroni_t_p_pCorr.tsv', 'w', newline='') as file:
    writer = csv.writer(file, delimiter='\t')
    for key, value in subfield_volumes_ipsVScontrols_ttestBonferroni_t_p_pCorr_dict.items():
        row = [key, value['t'], value['p'], value['pCorr_bonferroni']]
        writer.writerow(row)

