#Group displacement stats

import numpy as np
import pandas as pd
from scipy import stats #For t-test
import matplotlib.pyplot as plt
import os

#Inputs
###Displacement values LR tsv 
allSubjs_displacement_tsv = pd.read_csv(snakemake.params.group_LRcombined_CA1_tsv_subfield_displacement_spaceUnfold_mPD_atEachAP, sep='\t')
#  Remove [] from the displacement_data_list column and convert string to float
allSubjs_displacement_tsv['displacement_data_list'] = allSubjs_displacement_tsv['displacement_data_list'].str.strip('[]').apply(lambda x: float(x))

###FFT displacement values LR tsv
allSubjs_FFT_tsv = pd.read_csv(snakemake.params.group_LRcombined_CA1_tsv_subfield_FFT_displacement_spaceUnfold_mPD_atEachAP, sep='\t')
#  Remove [] from the FFT_displacement_data_list column and convert string to float
allSubjs_FFT_tsv['FFT_displacement_data_list'] = allSubjs_FFT_tsv['FFT_displacement_data_list'].str.strip('[]').apply(lambda x: float(x))

###PSD of FFT displacement values LR tsv
allSubjs_PSDofFFT_tsv = pd.read_csv(snakemake.params.group_LRcombined_CA1_tsv_subfield_FFT_PSD_displacement_spaceUnfold_mPD_atEachAP, sep='\t')
#  Remove [] from the PSDofFFT_displacement_data_list column and convert string to float
allSubjs_PSDofFFT_tsv['PSDofFFT_displacement_data_list'] = allSubjs_PSDofFFT_tsv['PSDofFFT_displacement_data_list'].str.strip('[]').apply(lambda x: float(x))

surface = snakemake.wildcards.surface

#Output_dirs
stats_allSubjs_displacement_tsv_dir = str(snakemake.output.stats_group_LRcombined_CA1_tsv_subfield_displacement_spaceUnfold_mPD_atEachAP_dir)
stats_allSubjs_FFT_displacement_tsv_dir = str(snakemake.output.stats_group_LRcombined_CA1_tsv_subfield_FFT_displacement_spaceUnfold_mPD_atEachAP_dir)
stats_allSubjs_PSDofFFT_displacement_tsv_dir = str(snakemake.output.stats_group_LRcombined_CA1_tsv_subfield_FFT_PSD_displacement_spaceUnfold_mPD_atEachAP_dir)


###1) Displacement values
#Separate grouped input into ipsilateral_TLE, contralateral_TLE, bilateral_TLE, and Controls
#For each of these groups, Group by participant_id to extract all the displacement_data_list values for each unique participant_id, so that a new df is created where each row contains a list of all the displacement_data_list values for that unique participant_id
#Ipsilateral (where tle_lateralization and hemi are the same)
ipsilateral_TLE_displacement_LR = allSubjs_displacement_tsv.query('tle_lateralization in ("R", "L") & hemi == tle_lateralization').copy()
ipsilateral_TLE_displacement_LR = ipsilateral_TLE_displacement_LR.groupby('participant_id')['displacement_data_list'].agg(list).reset_index()
print(ipsilateral_TLE_displacement_LR)
#Contralateral (where tle_lateralization is opposite of current hemi)
contralateral_TLE_displacement_LR = allSubjs_displacement_tsv.query('tle_lateralization in ("R", "L") & hemi != tle_lateralization').copy()
contralateral_TLE_displacement_LR = contralateral_TLE_displacement_LR.groupby('participant_id')['displacement_data_list'].agg(list).reset_index()
#Bilateral (have to group by participand_id AND hemi here, since both the L and R are used for each participant... hemi column is then dropped)
bilateral_TLE_displacement_LR = allSubjs_displacement_tsv.query('tle_lateralization == "BL"').copy()
bilateral_TLE_displacement_LR = bilateral_TLE_displacement_LR.groupby(['participant_id', 'hemi'])['displacement_data_list'].agg(list).reset_index().drop('hemi', axis=1)
#All affected hippocampi (ipsilateral_TLE + Bilateral_TLE)
allAffected_TLE_displacement_LR = pd.concat([ipsilateral_TLE_displacement_LR, bilateral_TLE_displacement_LR]).reset_index(drop=True)
#All TLE (ipsilateral_TLE_displacement + contralateral_TLE_displacement + bilateral_TLE_displacement)
all_TLE_displacement_LR = pd.concat([ipsilateral_TLE_displacement_LR, contralateral_TLE_displacement_LR, bilateral_TLE_displacement_LR]).reset_index(drop=True)
#Controls (have to group by participand_id AND hemi here, since both the L and R are used for each participant... hemi column is then dropped)
controls_displacement_LR = allSubjs_displacement_tsv.query('tle_lateralization == "Control"').copy()
controls_displacement_LR = controls_displacement_LR.groupby(['participant_id', 'hemi'])['displacement_data_list'].agg(list).reset_index().drop('hemi', axis=1)

#Calculate mean, std, sterror of displacement data across subjects for each group (230 values for each)
# Define a function to calculate mean, std, and stderr for each group
def calculate_group_stats(df, data_col_name):
    # Create an array from the displacement data list column
    data = np.array(df[data_col_name].tolist())
    # Calculate mean, std, and standard error of the mean
    mean = np.mean(data, axis=0)
    std = np.std(data, axis=0)
    stderr = std / np.sqrt(data.shape[0])
    # Return a dictionary of the calculated statistics
    return {'mean': mean, 'std': std, 'stderr': stderr}
ipsilateral_TLE_LR_stats_displacement = calculate_group_stats(ipsilateral_TLE_displacement_LR, 'displacement_data_list')
contralateral_TLE_LR_stats_displacement = calculate_group_stats(contralateral_TLE_displacement_LR, 'displacement_data_list')
bilateral_TLE_LR_stats_displacement = calculate_group_stats(bilateral_TLE_displacement_LR, 'displacement_data_list')
allAffected_TLE_LR_stats_displacement = calculate_group_stats(allAffected_TLE_displacement_LR, 'displacement_data_list')
all_TLE_LR_stats_displacement = calculate_group_stats(all_TLE_displacement_LR, 'displacement_data_list')
controls_LR_stats_displacement = calculate_group_stats(controls_displacement_LR, 'displacement_data_list')

###PLOTTING
#1) Plot means and std for each group's displacement_data_list (230 values for each group)
#2) Plot the std curves for each group's displacement_data_list
fig_displacement, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 12))
x_displacement = range(len(ipsilateral_TLE_LR_stats_displacement['mean']))
colors = ['red', 'blue', 'orange', 'brown', 'purple', 'green']
labels = ['Ipsilateral TLE', 'Contralateral TLE', 'Bilateral TLE', 'Ipsilateral + Bilateral TLE', 'Ipsilateral + Contralateral + Bilateral TLE', 'Controls']
for i, group in enumerate([ipsilateral_TLE_LR_stats_displacement, contralateral_TLE_LR_stats_displacement, bilateral_TLE_LR_stats_displacement, allAffected_TLE_LR_stats_displacement, all_TLE_LR_stats_displacement, controls_LR_stats_displacement]):
    ax1.errorbar(x_displacement, group['mean'], yerr=group['stderr'], label=labels[i], color=colors[i])
    ax2.plot(x_displacement, group['std'], label=labels[i], color=colors[i])
ax1.set_xticks([0, len(ipsilateral_TLE_LR_stats_displacement['mean'])])
ax1.set_xticklabels(['Anterior', 'Posterior'], fontsize=14)
ax1.set_xlabel("AP level", fontsize=14, labelpad=15)
ax1.set_ylabel("Group Mean Displacement", fontsize=14, labelpad=15)
ax1.set_title(str('Group Mean displacement at Each AP Level'), fontsize=16)
ax1.legend()
ax2.set_xticks([0, len(ipsilateral_TLE_LR_stats_displacement['mean'])])
ax2.set_xticklabels(['Anterior', 'Posterior'], fontsize=14)
ax2.set_xlabel("AP level", fontsize=14, labelpad=15)
ax2.set_ylabel("Group Mean Displacement", fontsize=14, labelpad=15)
ax2.set_title(str('Group Standard Deviation of Displacement at Each AP Level'), fontsize=16)
ax2.legend()
#plt.show()


###2) FFT of displacement values
#Separate grouped input into ipsilateral_TLE, contralateral_TLE, bilateral_TLE, and Controls
#For each of these groups, Group by participant_id to extract all the FFT_displacement_data_list values for each unique participant_id, so that a new df is created where each row contains a list of all the FFT_displacement_data_list values for that unique participant_id
#Ipsilateral (where tle_lateralization and hemi are the same)
ipsilateral_TLE_FFT_LR = allSubjs_FFT_tsv.query('tle_lateralization in ("R", "L") & hemi == tle_lateralization').copy() 
ipsilateral_TLE_FFT_LR = ipsilateral_TLE_FFT_LR.groupby('participant_id')['FFT_displacement_data_list'].agg(list).reset_index()
print(ipsilateral_TLE_FFT_LR)
#Contralateral (where tle_lateralization is opposite of current hemi)
contralateral_TLE_FFT_LR = allSubjs_FFT_tsv.query('tle_lateralization in ("R", "L") & hemi != tle_lateralization').copy()
contralateral_TLE_FFT_LR = contralateral_TLE_FFT_LR.groupby('participant_id')['FFT_displacement_data_list'].agg(list).reset_index()
#Bilateral (have to group by participand_id AND hemi here, since both the L and R are used for each participant)
bilateral_TLE_FFT_LR = allSubjs_FFT_tsv.query('tle_lateralization == "BL"').copy()
bilateral_TLE_FFT_LR = bilateral_TLE_FFT_LR.groupby(['participant_id', 'hemi'])['FFT_displacement_data_list'].agg(list).reset_index().drop('hemi', axis=1)
#All affected hippocampi (ipsilateral_TLE + Bilateral_TLE)
allAffected_TLE_FFT_LR = pd.concat([ipsilateral_TLE_FFT_LR, bilateral_TLE_FFT_LR]).reset_index(drop=True)
#All TLE (ipsilateral_TLE_FFT + contralateral_TLE_FFT + bilateral_TLE_FFT)
all_TLE_FFT_LR = pd.concat([ipsilateral_TLE_FFT_LR, contralateral_TLE_FFT_LR, bilateral_TLE_FFT_LR]).reset_index(drop=True)
#Controls (have to group by participand_id AND hemi here, since both the L and R are used for each participant)
controls_FFT_LR = allSubjs_FFT_tsv.query('tle_lateralization == "Control"').copy()
controls_FFT_LR = controls_FFT_LR.groupby(['participant_id', 'hemi'])['FFT_displacement_data_list'].agg(list).reset_index().drop('hemi', axis=1)

#Calculate mean, std, sterror of FFT data across subjects for each group (230 values for each)
ipsilateral_TLE_LR_stats_FFT = calculate_group_stats(ipsilateral_TLE_FFT_LR, 'FFT_displacement_data_list')
contralateral_TLE_LR_stats_FFT = calculate_group_stats(contralateral_TLE_FFT_LR, 'FFT_displacement_data_list')
bilateral_TLE_LR_stats_FFT = calculate_group_stats(bilateral_TLE_FFT_LR, 'FFT_displacement_data_list')
allAffected_TLE_LR_stats_FFT = calculate_group_stats(allAffected_TLE_FFT_LR, 'FFT_displacement_data_list')
all_TLE_LR_stats_FFT = calculate_group_stats(all_TLE_FFT_LR, 'FFT_displacement_data_list')
controls_LR_stats_FFT = calculate_group_stats(controls_FFT_LR, 'FFT_displacement_data_list')

# Define a function for absolute magnitude to calculate mean, std, and stderr for each group
def calculate_group_stats_abs_mag(df_FFT, stats_FFT):
    L = len(stats_FFT['mean'])
    T = 1.0/len(stats_FFT['mean'])
    # Create an array from the displacement data list column
    data = np.array(df_FFT[data_col_name].tolist())
    abs_mag = np.abs(data[:, :int(L/2)]) #Indexes each subject then takes np.abs of the first half of FFT data for each
    print(abs_mag[0], abs_mag.shape[0], abs_mag.shape[1])
    df_FFT['FFTabsmag_displacement_data_list'] = [list(np.ndarray.flatten(i)) for i in np.ndarray.reshape(abs_mag, (abs_mag.shape[0], abs_mag.shape[1]))]
    # Calculate mean, std, and standard error of the mean
    mean = np.mean(abs_mag, axis=0)
    std = np.std(abs_mag, axis=0)
    stderr = std / np.sqrt(abs_mag.shape[0]) 
    # Return the new df_FFT with absmag added and a dictionary of the calculated statistics
    return df_FFT, {'mean': mean, 'std': std, 'stderr': stderr}

ipsilateral_TLE_FFT_LR, ipsilateral_TLE_LR_stats_FFTmag = calculate_group_stats_abs_mag(ipsilateral_TLE_FFT_LR, 'FFT_displacement_data_list', ipsilateral_TLE_LR_stats_FFT)
contralateral_TLE_FFT_LR, contralateral_TLE_LR_stats_FFTmag = calculate_group_stats_abs_mag(contralateral_TLE_FFT_LR, 'FFT_displacement_data_list', contralateral_TLE_LR_stats_FFT)
bilateral_TLE_FFT_LR, bilateral_TLE_LR_stats_FFTmag = calculate_group_stats_abs_mag(bilateral_TLE_FFT_LR, 'FFT_displacement_data_list', bilateral_TLE_LR_stats_FFT)
allAffected_TLE_FFT_LR, allAffected_TLE_LR_stats_FFTmag = calculate_group_stats_abs_mag(allAffected_TLE_FFT_LR, 'FFT_displacement_data_list', allAffected_TLE_LR_stats_FFT)
all_TLE_FFT_LR, all_TLE_LR_stats_FFTmag = calculate_group_stats_abs_mag(all_TLE_FFT_LR, 'FFT_displacement_data_list', all_TLE_LR_stats_FFT)
controls_FFT_LR, controls_LR_stats_FFTmag = calculate_group_stats_abs_mag(controls_FFT_LR, 'FFT_displacement_data_list', controls_LR_stats_FFT)

###PLOTTING
#Plot FFT absolute magnitude means and std (Frequency Hz vs. Absolute Magnitude) for each group (115 values for each group)
L = L = len(ipsilateral_TLE_LR_stats_FFT['mean']) 
T = 1.0/len(ipsilateral_TLE_LR_stats_FFT['mean'])
xf = np.linspace(0.0, 1.0/(2.0*T), int(L/2))
mm_range_approx = 40
fig_FFT = plt.figure(figsize=(24,6))
colors = ['red', 'blue', 'orange', 'brown', 'purple', 'green']
labels = ['Ipsilateral TLE', 'Contralateral TLE', 'Bilateral TLE', 'Ipsilateral + Bilateral TLE', 'Ipsilateral + Contralateral + Bilateral TLE', 'Controls']
for i, group in enumerate([ipsilateral_TLE_LR_stats_FFTmag, contralateral_TLE_LR_stats_FFTmag, bilateral_TLE_LR_stats_FFTmag, allAffected_TLE_LR_stats_FFTmag, all_TLE_LR_stats_FFTmag, controls_LR_stats_FFTmag]):
    plt.errorbar(xf, 2.0/L * group['mean'], yerr=group['stderr'], label=labels[i], color=colors[i])
plt.xlabel('Frequency (Hz)')
plt.ylabel('Magnitude(abs)')
plt.title("Group Mean FFT Frequency (Hz) vs. Absolute Magnitude Displacement at Each AP Level", fontsize=16)
plt.legend()


###3) PSDofFFT of Displacement values
#Separate grouped input into ipsilateral_TLE, contralateral_TLE, bilateral_TLE, and Controls
#For each of these groups, Group by participant_id to extract all the PSDofFFT_displacement_data_list values for each unique participant_id, so that a new df is created where each row contains a list of all the PSDofFFT_displacement_data_list values for that unique participant_id
#Ipsilateral (where tle_lateralization and hemi are the same)
ipsilateral_TLE_PSDofFFT_LR = allSubjs_PSDofFFT_tsv.query('tle_lateralization in ("R", "L") & hemi == tle_lateralization').copy()
ipsilateral_TLE_PSDofFFT_LR = ipsilateral_TLE_PSDofFFT_LR.groupby('participant_id')['PSDofFFT_displacement_data_list'].agg(list).reset_index()
print(ipsilateral_TLE_PSDofFFT_LR)
#Contralateral (where tle_lateralization is opposite of current hemi)
contralateral_TLE_PSDofFFT_LR = allSubjs_PSDofFFT_tsv.query('tle_lateralization in ("R", "L") & hemi != tle_lateralization').copy()
contralateral_TLE_PSDofFFT_LR = contralateral_TLE_PSDofFFT_LR.groupby('participant_id')['PSDofFFT_displacement_data_list'].agg(list).reset_index()
#Bilateral (have to group by participand_id AND hemi here, since both the L and R are used for each participant)
bilateral_TLE_PSDofFFT_LR = allSubjs_PSDofFFT_tsv.query('tle_lateralization == "BL"').copy()
bilateral_TLE_PSDofFFT_LR = bilateral_TLE_PSDofFFT_LR.groupby(['participant_id', 'hemi'])['PSDofFFT_displacement_data_list'].agg(list).reset_index().drop('hemi', axis=1)
#All affected hippocampi (ipsilateral_TLE + Bilateral_TLE)
allAffected_TLE_PSDofFFT_LR = pd.concat([ipsilateral_TLE_PSDofFFT_LR, bilateral_TLE_PSDofFFT_LR]).reset_index(drop=True)
#All TLE (ipsilateral_TLE_PSDofFFT + contralateral_TLE_PSDofFFT + bilateral_TLE_PSDofFFT)
all_TLE_PSDofFFT_LR = pd.concat([ipsilateral_TLE_PSDofFFT_LR, contralateral_TLE_PSDofFFT_LR, bilateral_TLE_PSDofFFT_LR]).reset_index(drop=True)
#Controls (have to group by participand_id AND hemi here, since both the L and R are used for each participant)
controls_PSDofFFT_LR = allSubjs_PSDofFFT_tsv.query('tle_lateralization == "Control"').copy()
controls_PSDofFFT_LR = controls_PSDofFFT_LR.groupby(['participant_id', 'hemi'])['PSDofFFT_displacement_data_list'].agg(list).reset_index().drop('hemi', axis=1)

#Calculate mean, std, sterror of PSDofFFT data across subjects for each group (230 values for each)
ipsilateral_TLE_LR_stats_PSDofFFT = calculate_group_stats(ipsilateral_TLE_PSDofFFT_LR, 'PSDofFFT_displacement_data_list')
contralateral_TLE_LR_stats_PSDofFFT = calculate_group_stats(contralateral_TLE_PSDofFFT_LR, 'PSDofFFT_displacement_data_list')
bilateral_TLE_LR_stats_PSDofFFT = calculate_group_stats(bilateral_TLE_PSDofFFT_LR, 'PSDofFFT_displacement_data_list')
allAffected_TLE_LR_stats_PSDofFFT = calculate_group_stats(allAffected_TLE_PSDofFFT_LR, 'PSDofFFT_displacement_data_list')
all_TLE_LR_stats_PSDofFFT = calculate_group_stats(all_TLE_PSDofFFT_LR, 'PSDofFFT_displacement_data_list')
controls_LR_stats_PSDofFFT = calculate_group_stats(controls_PSDofFFT_LR, 'PSDofFFT_displacement_data_list')

#PLOTTING
#Plot PSD means and std for each group (Frequency Hz vs. Absolute PSD dB)
fig_PSDofFFT = plt.figure(figsize=(24,6))
colors = ['red', 'blue', 'orange', 'brown', 'purple', 'green']
labels = ['Ipsilateral TLE', 'Contralateral TLE', 'Bilateral TLE', 'Ipsilateral + Bilateral TLE', 'Ipsilateral + Contralateral + Bilateral TLE', 'Controls']
for i, group in enumerate([ipsilateral_TLE_LR_stats_PSDofFFT, contralateral_TLE_LR_stats_PSDofFFT, bilateral_TLE_LR_stats_PSDofFFT, allAffected_TLE_LR_stats_PSDofFFT, all_TLE_LR_stats_PSDofFFT, controls_LR_stats_PSDofFFT]):
    plt.errorbar(range(len(ipsilateral_TLE_LR_stats_PSDofFFT['mean'])), group['mean'], yerr=group['stderr'], label=labels[i], color=colors[i])
plt.xlabel('Frequency (Hz)')
plt.ylabel('PSD (dB)')
plt.title("Group Mean Power Spectrum Density of FFT Frequency (Hz) vs. Absolute PSD (dB) Displacement at Each AP Level", fontsize=16)
plt.legend()




#####Create output directories and save tsvs
#Groupnames
groupnames = ['ipsilateral_TLE', 'contralateral_TLE', 'bilateral_TLE', 'allAffected_TLE', 'all_TLE', 'controls']
#Displacement grouped data
displacement_dfs = [ipsilateral_TLE_displacement_LR, contralateral_TLE_displacement_LR, bilateral_TLE_displacement_LR, allAffected_TLE_displacement_LR, all_TLE_displacement_LR, controls_displacement_LR]
displacement_stats_dicts = [ipsilateral_TLE_LR_stats_displacement, contralateral_TLE_LR_stats_displacement, bilateral_TLE_LR_stats_displacement, allAffected_TLE_LR_stats_displacement, all_TLE_LR_stats_displacement, controls_LR_stats_displacement]
#FFT grouped data
FFT_dfs = [ipsilateral_TLE_FFT_LR, contralateral_TLE_FFT_LR, bilateral_TLE_FFT_LR, allAffected_TLE_FFT_LR, all_TLE_FFT_LR, controls_FFT_LR]
FFT_stats_dicts = [ipsilateral_TLE_LR_stats_FFT, contralateral_TLE_LR_stats_FFT, bilateral_TLE_LR_stats_FFT, allAffected_TLE_LR_stats_FFT, all_TLE_LR_stats_FFT, controls_LR_stats_FFT]
#PSDofFFT grouped data
PSDofFFT_dfs = [ipsilateral_TLE_PSDofFFT_LR, contralateral_TLE_PSDofFFT_LR, bilateral_TLE_PSDofFFT_LR, allAffected_TLE_PSDofFFT_LR, all_TLE_PSDofFFT_LR, controls_PSDofFFT_LR]
PSDofFFT_stats_dicts = [ipsilateral_TLE_LR_stats_PSDofFFT, contralateral_TLE_LR_stats_PSDofFFT, bilateral_TLE_LR_stats_PSDofFFT, allAffected_TLE_LR_stats_PSDofFFT, all_TLE_LR_stats_PSDofFFT, controls_LR_stats_PSDofFFT]


# Define the data types and their corresponding file names
#Grouped data df tsvs have columns participant_id and displacement_data_list(series) with rows=#participants
#  For FFT it has columns participant_id, FFT_displacement_data_list, and FFTabsmag_displacement_data_list
#  For PSDofFFT it has columns participant_id, PSDofFFT_displacement_data_list
#Gropued stats df tsvs have columns mean std sterror with rows=length displacement data (230)
def save_pd_tsvs(outdir, data_list_group, metric_name, stats=False):   #e.g. save_pd_tsvs(stats_allSubjs_displacement_tsv_dir, displacement_dfs, 'displacement_values')
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        for i,data in enumerate(data_list_group):
            if not stats:
                data.to_csv(outdir+'/group_surf-'+surface+'_'+groupnames[i]+'_'+metric_name+'.tsv', sep='\t', index=False)
            if stats:
                df_stats = pd.DataFrame.from_dict(data)
                df_stats.to_csv(outdir+'/group_surf-'+surface+'_'+groupnames[i]+'_'+metric_name+'.tsv', sep='\t', index=False)


#Displacement
save_pd_tsvs(stats_allSubjs_displacement_tsv_dir, displacement_dfs, 'displacement_values')
save_pd_tsvs(stats_allSubjs_displacement_tsv_dir, displacement_stats_dicts, 'mean_std_sterror_displacement', stats=True)
#FFT displacement
save_pd_tsvs(stats_allSubjs_FFT_displacement_tsv_dir, FFT_dfs, 'FFT_displacement_values')
save_pd_tsvs(stats_allSubjs_FFT_displacement_tsv_dir, FFT_stats_dicts, 'mean_std_sterror_FFT_displacement', stats=True)
#PSDofFFT displacement
save_pd_tsvs(stats_allSubjs_PSDofFFT_displacement_tsv_dir, PSDofFFT_dfs, 'PSDofFFT_displacement_values')
save_pd_tsvs(stats_allSubjs_PSDofFFT_displacement_tsv_dir, PSDofFFT_stats_dicts, 'mean_std_sterror_PSDofFFT_displacement', stats=True)

#Save figures
fig_displacement.savefig(stats_allSubjs_displacement_tsv_dir+'/allgroups_lineplot'+'surf-'+surface+'_TLE_MeanStd_and_std_line_displacement_values.png', dpi=300)
fig_FFT.savefig(stats_allSubjs_FFT_displacement_tsv_dir+'/allgroups_lineplot'+'surf-'+surface+'_TLE_MeanStd_FFT_frequency_vs_absMagnitude_displacement_values.png', dpi=300)
fig_PSDofFFT.savefig(stats_allSubjs_PSDofFFT_displacement_tsv_dir+'/allgroups_lineplot'+'surf-'+surface+'_TLE_MeanStd_PSDofFFT_frequency_vs_absPSD_displacement_values.png', dpi=300)
