# Additional hippunfold metrics to include in the pipeline. These are all analyzed in a similar way to displacement from the displacementT1w_to_displacementUnfold_nii rule onwards. These are only performed on midthickness surface as these are the current outputs of hippunfold

#Map metric to unfolded VOLUME space, using unfolded midsurf and constrained by unfolded inner/outer surfs
rule metricT1w_to_metricUnfold_nii:
    input:
        metric_T1w = config['shape_metric_midthickness_spaceT1w'], #Metric to map to vol
        gii_Unfolded_midthickness = expand(config['surf_spaceUnfold'], surface='midthickness', allow_missing=True), #Unfolded midsurf
        nii_Unfolded = config['refvol_spaceUnfold'] #Ref Unfolded vol to map metric to
    params:
        gii_Unfolded_inner = expand(config['surf_spaceUnfold'], surface='inner', allow_missing=True), #For -ribbon-constrained
        gii_Unfolded_outer = expand(config['surf_spaceUnfold'], surface='outer', allow_missing=True), #For -ribbon-constrained
    output:
        metric_Unfolded_midthickness_nii = 'work/hippunfold_morphometry/sub-{subject}/anat/sub-{subject}_hemi-{hemi}_space-unfolded_den-0p5mm_label-hipp_midthickness_{metric}.nii.gz'
    group: 'metric_subj_hemi'
    log: 'logs/hippunfold_morphometry/sub-{subject}/sub-{subject}_hemi-{hemi}_shape_{metric}_midthickness_T1w_xfm_spaceUnfold_nii.log'
    shell:
        'wb_command -metric-to-volume-mapping {input[0]} {input[1]} {input[2]} {output} -ribbon-constrained {params[0]} {params[1]} &> {log}'


#Metric analyses, outputs are as follows and all exist for each subfield and split into wholeHc(5% clipped from head and tail), HeadOnly, and BodyOnly (individual surfaces and hemis are processed in separate calls).
#Each AP's average (across PD) metric values (mean-centred) in each individual subfield (tsv and plot)
#fft of this (tsv and plot); plot excludes DC since data was mean-centred (saved in FFT outdir)
#psd of fft (tsv and plot); plot excludes DC since data was mean-centred (saved in FFT outdir)
rule metric_analyses_spaceUnfold:
    input:
        metric_Unfolded_midthickness_nii = rules.metricT1w_to_metricUnfold_nii.output,
        label_subfields_spaceUnfold_nii = expand(rules.subfieldsT1w_to_spaceUnfold_nii.output, surface='midthickness', allow_missing=True),
    params:
        subfield_names = subfields
    output:
        metric_subfields_analyses_Unfolded_to_smoothed_dir = directory('work/hippunfold_morphometry/sub-{subject}/{metric}_analyses/{hemi}_surf-midthickness/{metric}_mPD_at_eachAP/'),
        metric_subfields_analyses_Unfolded_to_smoothed_FFTandPSD_dir = directory('work/hippunfold_morphometry/sub-{subject}/{metric}_analyses/{hemi}_surf-midthickness/{metric}_mPD_at_eachAP_FFTandPSD/'),
    group: 'metric_subj_hemi'
    log: 'logs/hippunfold_morphometry/sub-{subject}/sub-{subject}_hemi-{hemi}_{metric}_analyses_spaceUnfold.log'
    script:
        '../scripts/metric_morphometry_midthickness_analysis.py'

#Group the tsv outputs of metric_analyses_spaceUnfold into combined tsvs to permit future analyses on the group data
#    The output tsvs are created for the following subgroups: ipsilateral_TLE, contralateral_TLE, bilateral_TLE, allAffected_TLE (ipsilateral + bilateral), all_TLE (ipsilateral, contralateral, bilateral), controls. 
#      Headers have 2 columns: participant_id and <name_of_metric> (which is a list with length 87 of the  <name_of_metric> data for that subject)
#        <name_of_metric> can be one of the following values 'mPDmetric_mAPcentered', 'FFT_mPDmetric_mAPcentered', or 'PSDofFFT_mPDmetric_mAPcentered'
#Currently only operating on CA1 and BodyOnly (see params)
rule group_metric_analyses_spaceUnfold:
    input:
        metric_subfields_analyses_Unfolded_to_smoothed_dir = expand(rules.metric_analyses_spaceUnfold.output[0], subject=subjects, hemi=hemis, allow_missing=True),
        metric_subfields_analyses_Unfolded_to_smoothed_FFTandPSD_dir = expand(rules.metric_analyses_spaceUnfold.output[1], subject=subjects, hemi=hemis, allow_missing=True),
    params:
        demographics_tsv = config['demographics_tsv'],
        #Actual inputs
        L_CA1_tsv_subfield_metric_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/{metric}_analyses/L_surf-midthickness/{metric}_mPD_at_eachAP/sub-{subject}_hemi-L_space-unfolded_subfield-CA1_BodyOnly_mPD{metric}_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
        L_CA1_tsv_FFT_subfield_metric_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/{metric}_analyses/L_surf-midthickness/{metric}_mPD_at_eachAP_FFTandPSD/sub-{subject}_hemi-L_space-unfolded_subfield-CA1_BodyOnly_FFT_mPD{metric}_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
        L_CA1_tsv_FFT_PSD_subfield_metric_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/{metric}_analyses/L_surf-midthickness/{metric}_mPD_at_eachAP_FFTandPSD/sub-{subject}_hemi-L_space-unfolded_subfield-CA1_BodyOnly_PSDofFFT_mPD{metric}_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
        R_CA1_tsv_subfield_metric_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/{metric}_analyses/R_surf-midthickness/{metric}_mPD_at_eachAP/sub-{subject}_hemi-R_space-unfolded_subfield-CA1_BodyOnly_mPD{metric}_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
        R_CA1_tsv_FFT_subfield_metric_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/{metric}_analyses/R_surf-midthickness/{metric}_mPD_at_eachAP_FFTandPSD/sub-{subject}_hemi-R_space-unfolded_subfield-CA1_BodyOnly_FFT_mPD{metric}_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
        R_CA1_tsv_FFT_PSD_subfield_metric_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/{metric}_analyses/R_surf-midthickness/{metric}_mPD_at_eachAP_FFTandPSD/sub-{subject}_hemi-R_space-unfolded_subfield-CA1_BodyOnly_PSDofFFT_mPD{metric}_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
    output:
        group_CA1_tsv_subfield_metric_spaceUnfold_mPD_atEachAP_dir = directory('results/hippunfold_morphometry/group_{metric}_analyses/surface-midthickness/'),
    threads: 8
    group: 'metric_groupSurfMid'
    log:
        'logs/hippunfold_morphometry/surf-midthickness_group_{metric}_analyses.log'
    script:
        '../scripts/group_metrics_morphometry_midthickness_analysis.py'


