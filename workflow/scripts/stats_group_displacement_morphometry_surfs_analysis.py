#Displacement morphometry statistical analyses:

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
import os


###Input dir
group_displacement_morphometry_groups_dir = str(snakemake.input)

###Params
#Subgroup names [ipsilateral_TLE, contralateral_TLE, bilateral_TLE, allAffected_TLE, all_TLE, controls]
group_names = snakemake.params[0]
#Displacement analysis types [mPDdisplacement_mAPcentered, FFT_mPDdisplacement_mAPcentered.tsv, PSDofFFT_mPDdisplacement_mAPcentered.tsv]
analysis_types = snakemake.params[1]

#Wildcards
surface = snakemake.wildcards.surface

#Output dir
stats_group_displacement_morphometry_groups_dir = str(snakemake.output)


#Define function to put all input subgroups for each analysis type into a dictionary of dfs
def create_groups_dict_dfs_displacement(groups_dir, names, analysis_type):
    groups_dict = {}
    for subgroup in names:
        f = groups_dir+'/group_surf-'+surface+'_'+subgroup+'_'+analysis_type+'.tsv'
        groups_dict[subgroup] = pd.read_csv(f, sep='\t')
    return groups_dict
#Create groups_dict_dfs_displacement into separate variables for each analysis_metric
groups_dict_dfs_mPDdisplacement = create_groups_dict_dfs_displacement(group_displacement_morphometry_groups_dir, group_names, str(analysis_types[0]))
groups_dict_dfs_FFT_mPDdisplacement = create_groups_dict_dfs_displacement(group_displacement_morphometry_groups_dir, group_names, str(analysis_types[1]))
groups_dict_dfs_PSDofFFT_mPDdisplacement = create_groups_dict_dfs_displacement(group_displacement_morphometry_groups_dir, group_names, str(analysis_types[2]))
print(groups_dict_dfs_mPDdisplacement)
print(groups_dict_dfs_FFT_mPDdisplacement)
print(groups_dict_dfs_PSDofFFT_mPDdisplacement)
#Define function to calculate energy of mPDdisplacement across bin_widths of 10%ofTotal to 50%ofTotal (for BodyOnly with total length of 87 this is 9 to 44 (rounded to int))
#For each group: Calculate the energy of bumpiness for each interval size by dividing the mPDdisplacement curve into intervals of using array_split, calculating the variance of each interval using var, summing up the variances using sum and dividing by the interval size to get the energy per interval. Energy per interval is appended to a list and are then summed up for all interval sizes to obtain the total energy of the bumpiness of the curve.
#These are placed into a dictionary for each subgroup where each row represents the total energy for that subject
def calc_group_energy_signal(groups_dict, min_bin_sz, max_bin_sz):
    groups_total_energy_dict = {}
    for subgroup, df in groups_dict.items():
        print(df['mPDdisplacement_mAPcentered'].values)
        subgroup_energies = []
        #Pull only the mPDdisplacement_values from the df, leaving a row for each subject in that group
        for interval_size in range(min_bin_sz-1, max_bin_sz):
            #For this group:  Split each subjects mPDdisplacement values into intervals where digitations are expected, use these intervals to calculate energy (sum of all variances across intervals) for each subject
            interval_split_mPDdisplacement = np.array_split(df['mPDdisplacement_mAPcentered'].values, len(df['mPDdisplacement_mAPcentered'].values) // interval_size)
            print(interval_split_mPDdisplacement.values)
            var_interval_split_mPDdisplacement = np.var(interval_split_mPDdisplacement, axis=1)
            energy_interval_split_mPDdisplacement = np.sum(var_interval_split_mPDdisplacement) / interval_size
            subgroup_energies.append(energy_interval_split_mPDdisplacement)
        groups_total_energy_dict[subgroup] = np.sum(subgroup_energies)
    return groups_total_energy_dict
groups_dict_total_energy_mPDdisplacement = calc_group_energy_signal(groups_dict_dfs_mPDdisplacement, 12, 29)
print(groups_dict_total_energy_mPDdisplacement)
#Create groups_total_energy_dict for mPDdisplacement signal
# Absolute max_bin_sz: Would never expect any digitations larger than 1/3 of the entire body (87/3 = 29datapoints)
# Absolute min_bin_sz: Would never expect to have more than 7 digitations in the body (87/7 = 12datapoints)
#    So lets test what group energy analyses look like with varying interval sizes between the following min_bin_sz, max_bin_sz; and plot them
#min_bin_szs = [12, 15, 18, 21]
#max_bin_szs = [29, 27, 25, 23]
#fig, axs = plt.subplots(len(min_bin_szs), len(max_bin_szs))
#for i, minbin in enumerate(min_bin_szs):
#    for j, maxbin in enumerate(min_bin_szs):
#        group_dict_total_energy_mPDdisplacement = calc_group_energy_signal(groups_dict_dfs_mPDdisplacement, minbin, maxbin)
#        print(group_dict_total_energy_mPDdisplacement)
#        mean_std_ipsilateral_TLE_total_energy_mPDdisplacement = {'mean': np.mean(group_dict_total_energy_mPDdisplacement['ipsilateral_TLE']), 'std': np.std(group_dict_total_energy_mPDdisplacement['ipsilateral_TLE'])}
#        print(mean_std_ipsilateral_TLE_total_energy_mPDdisplacement)
#        mean_std_contralateral_TLE_total_energy_mPDdisplacement = {'mean': np.mean(group_dict_total_energy_mPDdisplacement['contralateral_TLE']), 'std': np.std(group_dict_total_energy_mPDdisplacement['contralateral_TLE'])}
#        print(mean_std_contralateral_TLE_total_energy_mPDdisplacement)
#        mean_std_controls_total_energy_mPDdisplacement = {'mean': np.mean(group_dict_total_energy_mPDdisplacement['controls']), 'std': np.std(group_dict_total_energy_mPDdisplacement['controls'])}
#        print(mean_std_controls_total_energy_mPDdisplacement)
#        x = np.arange(3) #3 Groups ipsilateral_TLE, contralateral_TLE, controls
#        axs[i][j].bar(x, [mean_std_ipsilateral_TLE_total_energy_mPDdisplacement['mean'], mean_std_contralateral_TLE_total_energy_mPDdisplacement['mean'], mean_std_controls_total_energy_mPDdisplacement['mean']], yerr=[mean_std_ipsilateral_TLE_total_energy_mPDdisplacement['std'], mean_std_contralateral_TLE_total_energy_mPDdisplacement['std'], mean_std_controls_total_energy_mPDdisplacement['std']], align='center', alpha=0.5)
#        axs[i][j].set_xticks(x)
#        axs[i][j].set_xtick_labels(['ipsilateral_TLE', 'contralateral_TLE', 'controls'])
#        axs[i][j].set_ylabel('Group mean total energy mPD_displacement')
#        axs[i][j].set_title('Group mean/std total energy mPD_displacement with interval_range ['+minbin+' to '+maxbin+']')
#        
##Save figure (each subplot is for a different interval_range) for each subgroups' mean_total_energy_mPDdisplacement across all subjects for that group
#if not os.path.exists(stats_group_displacement_morphometry_groups_dir):
#    os.mkdir(stats_group_displacement_morphometry_groups_dir)
#
#fig.savefig(stats_group_displacement_morphometry_groups_dir+'/groups_meanStd_totalEnergy_mPDdisplacement_allTESTintervalRanges.png')





