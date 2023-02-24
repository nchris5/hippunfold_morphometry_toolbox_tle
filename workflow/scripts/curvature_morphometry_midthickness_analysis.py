#Curvature

import numpy as np
import pandas as pd
import nibabel as nib
import statistics
import os
import copy
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

#Inputs
#1) Curvature from unsmoothed to smoothed surface in unfolded nii space
curvature_Unfolded_to_smoothed_nii = snakemake.input[0]
#2) Subfield labels in unfolded nii space)
labelSubfields_Unfolded_nii = snakemake.input[1]

#Params:
#1) subfield names (Subiculum, CA1, CA2, CA3, CA4)
subfield_names = snakemake.params.subfield_names
#2) Methods to extract different AP lengths from the whole hippocampus
extraction_method = ['remove5pHeadTail', 'HeadOnly', 'BodyOnly']
#3) Surface used
surface = 'midthickness'

#Wildcards: subject and hemispehere
subj = str(snakemake.wildcards.subject)
hemi = str(snakemake.wildcards.hemi)
print('Subject =', subj)
print('Hemi =', hemi)

#Outputs:
outdir_curvature_mPD_at_eachAP_data = str(snakemake.output[0])
outdir_curvature_mPD_at_eachAP_FFTandPSD = str(snakemake.output[1])


#Loop all the following operations for each subfield (Subiculum=1, CA1=2, CA2=3, CA3=4, CA4=5))
for sf_idx,sf_name in enumerate(subfield_names):
    #Load data from inputs (curvature_nii and sfLabels_nii)
    curvature_nii = nib.load(curvature_Unfolded_to_smoothed_nii).get_fdata()
    sfLabels_nii = nib.load(labelSubfields_Unfolded_nii).get_fdata()

    #Only need a single AP*PD slice across IO, so we take the midthickness (only the 1st and 16th (inner and outer) slices are 0, the rest are equal)
    curvature_nii_midIO = curvature_nii[:,:,7]
    sfLabels_midIO = sfLabels_nii[:,:,7]

    #Get subfield indices in unfolded volume space at the midIO (in actuality the notSubfield indices will be used to set anything not in current subfield to NaN)
    sfLabels_midIO_notIdx_sf = np.where(sfLabels_midIO != sf_idx+1)
    #Use these notSubield indices to set the curvature_nii at the midIO to NaN outside of current subfield
    curvature_nii_midIO_sf = curvature_nii_midIO #Initialize new variable
    curvature_nii_midIO_sf[sfLabels_midIO_notIdx_sf] = np.nan
    
    #Loop through AP curvature values in current subfield and take the average across PD at each level of AP
    curvature_nii_midIO_sf_mPD_eachAP = [np.nanmean(AP) for AP in curvature_nii_midIO_sf]

    #Loop through 'remove5pHeadTail', 'HeadOnly', 'BodyOnly' and do the remainder of analyses for each of these individually
    #Also apply np.abs at this stage since one hemi is negative
    for method in extraction_method:
        if method == 'remove5pHeadTail':
            curvature_nii_midIO_sf_mPD_extractAP = np.abs(curvature_nii_midIO_sf_mPD_eachAP[12:-13])
            #Center the curvature values around the mean curvature
            curvature_nii_midIO_sf_mPD_extractAP_mAPcentered = curvature_nii_midIO_sf_mPD_extractAP - np.mean(curvature_nii_midIO_sf_mPD_extractAP)
        elif method == 'HeadOnly':
            curvature_nii_midIO_sf_mPD_extractAP = np.abs(curvature_nii_midIO_sf_mPD_eachAP[12:int(round(len(curvature_nii_midIO_sf_mPD_eachAP)/3))])
            curvature_nii_midIO_sf_mPD_extractAP_mAPcentered = curvature_nii_midIO_sf_mPD_extractAP - np.mean(curvature_nii_midIO_sf_mPD_extractAP)
        elif method == 'BodyOnly':
            curvature_nii_midIO_sf_mPD_extractAP = np.abs(curvature_nii_midIO_sf_mPD_eachAP[int(round(len(curvature_nii_midIO_sf_mPD_eachAP)/3)-1):int(round(len(curvature_nii_midIO_sf_mPD_eachAP)*2/3))])
            curvature_nii_midIO_sf_mPD_extractAP_mAPcentered = curvature_nii_midIO_sf_mPD_extractAP - np.mean(curvature_nii_midIO_sf_mPD_extractAP)

        #Fourier analysis (FFT) on Mean PD absolute and mean centered curvature at each level of AP
        fft_curvature_nii_midIO_sf_mPD_extractAP_mAPcentered = np.fft.fft(curvature_nii_midIO_sf_mPD_extractAP_mAPcentered)

        L = len(curvature_nii_midIO_sf_mPD_extractAP)
        T = 1.0/len(curvature_nii_midIO_sf_mPD_extractAP)
        #t = np.linspace(0.0, T*L, L) #Not used, but if different L and T are used, can create an x-axis for arbitrary L length 
        #xf = np.linspace(0.0, 1.0/(2.0*T), int(L/2)) #Not used, but this is used to plot arbitrary frequency instead of real world mm values

        #Get frequencies in mm (mm_range_approx is 40 * ratio of original length (256) to clipped length (e.g. 87 for BodyOnly)
        mm_range_approx = 40*(len(curvature_nii_midIO_sf_mPD_extractAP)/len(curvature_nii_midIO_sf_mPD_eachAP))
        freq_mm = np.fft.fftfreq(L, d=1/(1/(mm_range_approx/L)))
        #Get Absolute value of freq_mm for plotting (this becomes x-axis)
        xf_freq_mm_abs = np.abs(freq_mm[:len(freq_mm)//2])
        #Convert freq_mm to the real length in mm that each freq represents (1/freq_mm), replace inf (0) with 'DC', and round to 3 decimal places.
        #This becomes x-axis labels)
        with np.errstate(divide='ignore', invalid='ignore'): #Hide the divide by zero outerr printing for freq=0 (since this crowds the output of snakemake). This is replaced with the string 'DC' later on, and not included when plotting FFTmag and PSD since DC has been mean centred
            xf_freqlabels_mm_real = np.abs(1/freq_mm[:len(freq_mm)//2]).round(3).astype(str)
        xf_freqlabels_mm_real[0] = 'DC'

        #Get absolute magnitudes
        magnitude = np.abs(fft_curvature_nii_midIO_sf_mPD_extractAP_mAPcentered[:int(L/2)])

        #NOT USED: Identify the dominant frequencies (peaks)
        peaks, _ = find_peaks(magnitude[:L//2])
        frequencies = peaks * T / L

        #Power Spectrum Density (PSD), measure of power persent in signal as a function of frequency (mag**2 of FFT)
        psd = np.abs((fft_curvature_nii_midIO_sf_mPD_extractAP_mAPcentered[:int(L/2)])**2)
        #Convert PSD to dB
        psd = 10*np.log10(psd)



        #Plotting
        #1) Plot Mean curvature across PD at each level of AP
        fig_curvature = plt.figure(figsize=(24,6))
        x = np.array(range(len(curvature_nii_midIO_sf_mPD_extractAP_mAPcentered)))
        plt.plot(x, curvature_nii_midIO_sf_mPD_extractAP_mAPcentered, '-o', markersize=2, label=str(hemi+' hemi'), color='blue')
        plt.xlim(0, len(curvature_nii_midIO_sf_mPD_extractAP_mAPcentered))
        plt.xticks([0, len(curvature_nii_midIO_sf_mPD_extractAP_mAPcentered)], ['Anterior', 'Posterior'], fontsize=14)
        plt.xlabel("AP level", fontsize=14, labelpad=15)
        plt.ylabel("Mean PD Curvature to Smoothed Surface", fontsize=14, labelpad=15)
        plt.title(str('Sub-'+subj+' '+hemi+' Hemisphere '+surface+' surface '+sf_name+' '+method+' Mean PD Curvature-to-Smoothed-Surface at Each AP Level'), fontsize=16)
        plt.legend()


        #2a) Plot FFT curvature absolute magnitude mPD at each level of AP (rescaled to mm)
        #Note that x-axis labels are the real world mm length that these peaks correspond to (1/freq_mm)
        #When plotting we skip x=0 since the DC (mean is 0)
        fig_FFTmag_and_PSD, (ax_FFTmag, ax_PSD) = plt.subplots(2, 1, figsize=(16,12))
        ax_FFTmag.plot(np.arange(1,len(xf_freq_mm_abs)), magnitude[1:], label=str(hemi+' hemi'), color='blue')
        ax_FFTmag.set_xticks(np.arange(1,len(xf_freq_mm_abs))) #Adds every value as an xtick
        ax_FFTmag.set_xticklabels(xf_freqlabels_mm_real[1:], rotation='vertical')
        ax_FFTmag.set_ylabel('Magnitude (abs)')
        ax_FFTmag.set_title(str('Sub-'+subj+' '+hemi+' Hemisphere '+surface+' surface '+sf_name+' '+method+' FFT Corresponding Length (mm) of Frequencies vs. Absolute Magnitude Curvature at Each AP Level'), fontsize=16)
        ax_FFTmag.legend()
        
        #2b) Plot PSD curvature absolute (dB)
        ax_PSD.plot(np.arange(1,len(xf_freq_mm_abs)), psd[1:], label=str(hemi+' hemi'), color='orange')
        ax_PSD.set_xlabel('Length (mm) of each frequency component (1/freq_mm)')
        ax_PSD.set_xticks(np.arange(1,len(xf_freq_mm_abs))) #Adds every value as an x-tick
        ax_PSD.set_xticklabels(xf_freqlabels_mm_real[1:], rotation='vertical')
        ax_PSD.set_ylabel('PSD (dB)')
        ax_PSD.set_title("Power Spectrum Density of FFT Corresponding Length (mm) of Frequencies vs. Absolute PSD (dB) Curvature at Each AP Level", fontsize=16)
        ax_PSD.legend()


        #Save output tsvs and figures
        if not os.path.exists(outdir_curvature_mPD_at_eachAP_data):
            os.mkdir(outdir_curvature_mPD_at_eachAP_data)
        if not os.path.exists(outdir_curvature_mPD_at_eachAP_FFTandPSD):
            os.mkdir(outdir_curvature_mPD_at_eachAP_FFTandPSD)

        #1) Save each AP's average (across PD) absolute and mean centered curvature values (also save non-mean centered) in current subfield (with current extraction method)
        np.savetxt(outdir_curvature_mPD_at_eachAP_data+'/sub-'+subj+'_hemi-'+hemi+'_space-unfolded_subfield-'+sf_name+'_'+method+'_mPDcurvature_at_each_AP.tsv', curvature_nii_midIO_sf_mPD_extractAP, delimiter='\t')
        np.savetxt(outdir_curvature_mPD_at_eachAP_data+'/sub-'+subj+'_hemi-'+hemi+'_space-unfolded_subfield-'+sf_name+'_'+method+'_mPDcurvature_mAPcentered_at_each_AP.tsv', curvature_nii_midIO_sf_mPD_extractAP_mAPcentered, delimiter='\t')
        #2) Save FFT and PSD of each AP's average (across PD) curvature values in current subfield
        np.savetxt(outdir_curvature_mPD_at_eachAP_FFTandPSD+'/sub-'+subj+'_hemi-'+hemi+'_space-unfolded_subfield-'+sf_name+'_'+method+'_FFT_mPDcurvature_mAPcentered_at_each_AP.tsv', fft_curvature_nii_midIO_sf_mPD_extractAP_mAPcentered, delimiter='\t')
        np.savetxt(outdir_curvature_mPD_at_eachAP_FFTandPSD+'/sub-'+subj+'_hemi-'+hemi+'_space-unfolded_subfield-'+sf_name+'_'+method+'_FFTfreqMM_mPDcurvature_mAPcentered_at_each_AP.tsv', freq_mm, delimiter='\t')
        np.savetxt(outdir_curvature_mPD_at_eachAP_FFTandPSD+'/sub-'+subj+'_hemi-'+hemi+'_space-unfolded_subfield-'+sf_name+'_'+method+'_PSDofFFT_mPDcurvature_mAPcentered_at_each_AP.tsv', psd, delimiter='\t')
        #3) Save and close figures
        fig_curvature.savefig(outdir_curvature_mPD_at_eachAP_data+'/sub-'+subj+'_hemi-'+hemi+'_subfield-'+sf_name+'_'+method+'_mPDcurvature_at_each_AP.png')
        fig_FFTmag_and_PSD.savefig(outdir_curvature_mPD_at_eachAP_FFTandPSD+'/sub-'+subj+'_hemi-'+hemi+'_subfield-'+sf_name+'_'+method+'_FFTmag_and_PSD_mPDcurvature_at_each_AP.png')
        plt.close(fig_curvature)
        plt.close(fig_FFTmag_and_PSD)