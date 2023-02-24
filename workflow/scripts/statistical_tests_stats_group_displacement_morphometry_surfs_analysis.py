#Statistical tests group displacement stats
import numpy as np
import pandas as pd
from scipy import stats #For any t-test
import matplotlib.pyplot as plt
import os
from scipy.stats import levene #Equality of variance (p < 0.05 means unequal variances)
from scipy.stats import mannwhitneyu #Non-equal-variance t-test
from scipy.signal import welch #For extracting specific frequency bands
from scipy.stats import f_oneway #Non-equal variance ANOVA (welch)
from statsmodels.stats.multicomp import pairwise_tukeyhsd #For post-hoc ANOVA analysis


groupnames = ['ipsilateral_TLE', 'contralateral_TLE', 'allAffected_TLE', 'controls']
surface = snakemake.wildcards.surface
#Inputs
indir_displacement = str(snakemake.input.stats_group_LRcombined_CA1_tsv_subfield_displacement_spaceUnfold_mPD_atEachAP_dir)
indir_FFT = str(snakemake.input.stats_group_LRcombined_CA1_tsv_subfield_FFT_displacement_spaceUnfold_mPD_atEachAP_dir)
indir_PSDofFFT = str(snakemake.input.stats_group_LRcombined_CA1_tsv_subfield_FFT_PSD_displacement_spaceUnfold_mPD_atEachAP_dir)

#1) Grouped data df tsvs have columns participant_id and displacement_data_list(series) with rows=#participants
#Displacement
ipsilateral_TLE_displacement = pd.read_csv(indir_displacement+'/group_surf-'+surface+'_ipsilateral_TLE_displacement_values.tsv', sep='\t')
contralateral_TLE_displacement = pd.read_csv(indir_displacement+'/group_surf-'+surface+'_contralateral_TLE_displacement_values.tsv', sep='\t')
allAffected_TLE_displacement = pd.read_csv(indir_displacement+'/group_surf-'+surface+'_allAffected_TLE_displacement_values.tsv', sep='\t') 
controls_displacement = pd.read_csv(indir_displacement+'/group_surf-'+surface+'_controls_displacement_values.tsv', sep='\t')
#FFT
ipsilateral_TLE_FFT = pd.read_csv(indir_FFT+'/group_surf-'+surface+'_ipsilateral_TLE_FFT_displacement_values.tsv', sep='\t')
contralateral_TLE_FFT = pd.read_csv(indir_FFT+'/group_surf-'+surface+'_contralateral_TLE_FFT_displacement_values.tsv', sep='\t')
allAffected_TLE_FFT = pd.read_csv(indir_FFT+'/group_surf-'+surface+'_allAffected_TLE_FFT_displacement_values.tsv', sep='\t')
controls_FFT = pd.read_csv(indir_FFT+'/group_surf-'+surface+'_controls_FFT_displacement_values.tsv', sep='\t')
#PSDofFFT
ipsilateral_TLE_PSDofFFT = pd.read_csv(indir_PSDofFFT+'/group_surf-'+surface+'_ipsilateral_TLE_PSDofFFT_displacement_values.tsv', sep='\t')
contralateral_TLE_PSDofFFT = pd.read_csv(indir_PSDofFFT+'/group_surf-'+surface+'_contralateral_TLE_PSDofFFT_displacement_values.tsv', sep='\t')
allAffected_TLE_PSDofFFT = pd.read_csv(indir_PSDofFFT+'/group_surf-'+surface+'_allAffected_TLE_PSDofFFT_displacement_values.tsv', sep='\t')
controls_PSDofFFT = pd.read_csv(indir_PSDofFFT+'/group_surf-'+surface+'_controls_PSDofFFT_displacement_values.tsv', sep='\t')

#2) Gropued stats df tsvs have columns mean std sterror with rows=length displacement data (230)
#Displacement
ipsilateral_TLE_displacement_stats = pd.read_csv(indir_displacement+'/group_surf-'+surface+'_ipsilateral_TLE_mean_std_sterror_displacement', sep='\t')
contralateral_TLE_displacement_stats = pd.read_csv(indir_displacement+'/group_surf-'+surface+'_contralateral_TLE_mean_std_sterror_displacement', sep='\t')
allAffected_TLE_displacement_stats = pd.read_csv(indir_displacement+'/group_surf-'+surface+'_allAffected_TLE_mean_std_sterror_displacement', sep='\t')
controls_displacement_stats = pd.read_csv(indir_displacement+'/group_surf-'+surface+'_controls_mean_std_sterror_displacement', sep='\t')
#FFT
ipsilateral_TLE_FFT_stats = pd.read_csv(indir_FFT+'/group_surf-'+surface+'_ipsilateral_TLE_mean_std_sterror_FFT_displacement', sep='\t')
contralateral_TLE_FFT_stats = pd.read_csv(indir_FFT+'/group_surf-'+surface+'_contralateral_TLE_mean_std_sterror_FFT_displacement', sep='\t')
allAffected_TLE_FFT_stats = pd.read_csv(indir_FFT+'/group_surf-'+surface+'_allAffected_TLE_mean_std_sterror_FFT_displacement', sep='\t')
controls_FFT_stats = pd.read_csv(indir_FFT+'/group_surf-'+surface+'_controls_mean_std_sterror_FFT_displacement', sep='\t')
#PSDofFFT
ipsilateral_TLE_PSDofFFT_stats = pd.read_csv(indir_FFT+'/group_surf-'+surface+'_ipsilateral_TLE_mean_std_sterror_FFT_PSDofFFTdisplacement', sep='\t')
contralateral_TLE_PSDofFFT_stats = pd.read_csv(indir_FFT+'/group_surf-'+surface+'_contralateral_TLE_mean_std_sterror_PSDofFFT_displacement', sep='\t')
allAffected_TLE_PSDofFFT_stats = pd.read_csv(indir_FFT+'/group_surf-'+surface+'_allAffected_TLE_mean_std_sterror_PSDofFFT_displacement', sep='\t')
controls_PSDofFFT_stats = pd.read_csv(indir_FFT+'/group_surf-'+surface+'_controls_mean_std_sterror_PSDofFFT_displacement', sep='\t')


#Outputs 
#
#


###1) Displacement Data Levene for equality of variances
# Function to compute the variances of each group
def calculate_group_var(df, data_col_name):
    # Create an array from the displacement data list column
    data = np.array(df[data_col_name].tolist())
    # Calculate mean, std, and standard error of the mean
    var = np.var(data, axis=0)
    # Return a dictionary of the calculated statistics
    return var
ipsilateral_TLE_displacement_stats['var'] = calculate_group_var(ipsilateral_TLE_displacement_stats, 'displacement_data_list')
contralateral_TLE_displacement_stats['var'] = calculate_group_var(contralateral_TLE_displacement_stats, 'displacement_data_list')
allAffected_TLE_displacement_stats['var'] = calculate_group_var(allAffected_TLE_displacement_stats, 'displacement_data_list')
controls_displacement_stats['var'] = calculate_group_var(controls_displacement_stats, 'displacement_data_list')
#Perform Levene's test ipsilateral vs. contralateral
ipsilateral_vs_contralateral_displacement_Levene_statistic, ipsilateral_vs_contralateral_displacement_Levene_p_value = levene(ipsilateral_TLE_displacement_stats['var'], contralateral_TLE_displacement_stats['var'])
print(f"Levene Ipsilateral vs Contralateral statistic: {ipsilateral_vs_contralateral_displacement_Levene_statistic:.3f}", f"p-value: {ipsilateral_vs_contralateral_displacement_Levene_p_value.3f}")
#Ipsilateral vs. Control
ipsilateral_vs_controls_displacement_Levene_statistic, ipsilateral_vs_controls_displacement_Levene_p_value = levene(ipsilateral_TLE_displacement_stats['var'], controls_displacement_stats['var'])
print(f"Levene Ipsilateral vs Controls statistic: {ipsilateral_vs_controls_displacement_Levene_statistic:.3f}", f"p-value: {ipsilateral_vs_controls_displacement_Levene_p_value.3f}")


###2) Displacement Data Mann-Whitney-U
ipsilateral_vs_contralateral_displacement_MWU_statistic, ipsilateral_vs_contralateral_displacement_MWU_p_value = mannwhitneyu(ipsilateral_TLE_displacement['displacement_data_list'], contralateral_TLE_displacement['displacement_data_list'])
print(f"Mann-Whitney-U Ipsilateral vs Contralateral statistic: {ipsilateral_vs_contralateral_displacement_MWU_statistic:.3f}", f"p-value: {ipsilateral_vs_contralateral_displacement_MWU_p_value.3f}")
#Ipsilateral vs. Control
ipsilateral_vs_controls_displacement_MWU_statistic, ipsilateral_vs_controls_displacement_MWU_p_value = mannwhitneyu(ipsilateral_TLE_displacement['displacement_data_list'], controls_displacement['displacement_data_list'])
print(f"Mann-Whitney-U Ipsilateral vs Controls statistic: {ipsilateral_vs_controls_displacement_MWU_statistic:.3f}", f"p-value: {ipsilateral_vs_controls_displacement_MWU_p_value.3f}")


###3) Displacement Data std t-test to compare the distributions of standard deviation values
ipsilateral_vs_contralateral_displacement_ttestStd_statistic, ipsilateral_vs_contralateral_displacement_ttestStd_p_value = stats.ttest_ind(ipsilateral_TLE_displacement_stats['std'], contralateral_TLE_displacement_stats['std'])
print(f"T-test on standard-deviation Ipsilateral vs Contralateral statistic: {ipsilateral_vs_contralateral_displacement_ttestStd_statistic:.3f}", f"p-value: {ipsilateral_vs_contralateral_displacement_ttestStd_p_value.3f}")
#Ipsilateral vs. Control
ipsilateral_vs_controls_displacement_ttestStd_statistic, ipsilateral_vs_controls_displacement_ttestStd_p_value = stats.ttest_ind(ipsilateral_TLE_displacement_stats['std'], controls_displacement_stats['std'])
print(f"T-test on standard-deviation Ipsilateral vs Controls statistic: {ipsilateral_vs_controls_displacement_ttestStd_statistic:.3f}", f"p-value: {ipsilateral_vs_controls_displacement_ttestStd_p_value.3f}")


###4) Welch method PSD calculation filtered for delta band (lower frequency variations 0.5-4) and subsequent Welch ANOVA (unequal variance) between ipsilateral vs. contralateral vs. controls
ipsilateral_TLE_welchPSD_delta = np.zeros((ipsilateral_TLE_displacement.shape[0], len(ipsilateral_TLE_displacement_stats['mean']/2))
for i in range(ipsilateral_TLE_displacement.shape[0]):
    ipsilateral_freq, ipsilateral_welchPSD = welch(ipsilateral_TLE_displacement['displacement_data_list'][i], nperseg=128)
    ipsilateral_TLE_welchPSD_delta[i] = ipsilateral_welchPSD[(ipsilateral_freq, >= 0.5) & (ipsilateral_freq <= 4)]

contralateral_TLE_welchPSD_delta = np.zeros((contralateral_TLE_displacement.shape[0], len(contralateral_TLE_displacement_stats['mean']/2))
for i in range(contralateral_TLE_displacement.shape[0]):
    contralateral_freq, contralateral_welchPSD = welch(contralateral_TLE_displacement['displacement_data_list'][i], nperseg=128)
    contralateral_TLE_welchPSD_delta[i] = contralateral_welchPSD[(contralateral_freq, >= 0.5) & (contralateral_freq <= 4)]

controls_welchPSD_delta = np.zeros((controls_displacement.shape[0], len(controls_displacement_stats['mean']/2))
for i in range(controls_displacement.shape[0]):
    controls_freq, controls_welchPSD = welch(controls_displacement['displacement_data_list'][i], nperseg=128)
    controls_welchPSD_delta[i] = controls_welchPSD[(controls_freq, >= 0.5) & (controls_freq <= 4)]

# Select PSD values for delta band
delta_ipsilateral_TLE_welchPSD_delta = ipsilateral_TLE_welchPSD_delta[:, :]
delta_contralateral_TLE_welchPSD_delta = contralateral_TLE_welchPSD_delta[:, :]
delta_controls_welchPSD_delta = controls_welchPSD_delta[:, :]

# Compare Welch delta PSD values between groups ANOVA
ipsilateral_vs_contralateral_vs_contols_welchfANOVA_statistic, ipsilateral_vs_contralateral_vs_contols_welchfANOVA_p_value = f_oneway(delta_ipsilateral_TLE_welchPSD_delta, delta_contralateral_TLE_welchPSD_delta, delta_controls_welchPSD_delta)
print(f"ipsilateral vs contralateral vs controls welch f_oneway ANOVA p-value for delta band: {ipsilateral_vs_contralateral_vs_contols_welchfANOVA_statistic:.3f}", f"p-value: {ipsilateral_vs_contralateral_vs_contols_welchfANOVA_p_value.3f}")

###4b) Post-hoc tukey-HSD of welch delta ANOVA
#Concat groups into a single list
delta_ipsilateral_contralateral_controls_delta_concat = np.concatenate([delta_ipsilateral_TLE_welchPSD_delta, delta_contralateral_TLE_welchPSD_delta, delta_controls_welchPSD_delta])
# Create group labels
delta_ipsilateral_contralateral_controls_delta_labels = ['ipsilateral_TLE'] * len(delta_ipsilateral_TLE_welchPSD_delta) + ['contralateral_TLE'] * len(delta_contralateral_TLE_welchPSD_delta) + ['controls'] * len(delta_controls_welchPSD_delta)

# Perform Tukey's HSD test
delta_ipsilateral_contralateral_controls_delta_tukeyHSD_results = pairwise_tukeyhsd(delta_ipsilateral_contralateral_controls_delta_concat, delta_ipsilateral_contralateral_controls_delta_labels, alpha=0.05)
print('Tukey HSD delta_ipsilateral_contralateral_controls_delta_tukeyHSD_results: ', delta_ipsilateral_contralateral_controls_delta_tukeyHSD_results)



