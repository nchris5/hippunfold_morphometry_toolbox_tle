#
#Map curvature to unfolded VOLUME space, using unfolded midsurf and constrained by unfolded inner/outer surfs
rule curvatureT1w_to_curvatureUnfold_nii:
    input:
        curvature_T1w = config['shape_curvature_midthickness_spaceT1w'], #Metric to map to vol
        gii_Unfolded_midthickness = expand(config['surf_spaceUnfold'], surface='midthickness', allow_missing=True), #Unfolded midsurf
        nii_Unfolded = config['refvol_spaceUnfold'] #Ref Unfolded vol to map metric to
    params:
        gii_Unfolded_inner = expand(config['surf_spaceUnfold'], surface='inner', allow_missing=True), #For -ribbon-constrained
        gii_Unfolded_outer = expand(config['surf_spaceUnfold'], surface='outer', allow_missing=True), #For -ribbon-constrained
    output:
        curvature_Unfolded_midthickness_nii = 'work/hippunfold_morphometry/sub-{subject}/anat/sub-{subject}_hemi-{hemi}_space-unfolded_den-0p5mm_label-hipp_midthickness_curvature.nii.gz'
    group: 'subjLR'
    log: 'logs/hippunfold_morphometry/sub-{subject}/sub-{subject}_hemi-{hemi}_shape_curvature_midthickness_T1w_xfm_spaceUnfold_nii.log'
    shell:
        'wb_command -metric-to-volume-mapping {input[0]} {input[1]} {input[2]} {output} -ribbon-constrained {params[0]} {params[1]} &> {log}'


#Curvature analyses, outputs are as follows and all exist for each subfield and split into wholeHc(5% clipped from head and tail), HeadOnly, and BodyOnly (individual surfaces and hemis are processed in separate calls).
#Each AP's average (across PD) curvature values (mean-centred) in each individual subfield (tsv and plot)
#fft of this (tsv and plot); plot excludes DC since data was mean-centred (saved in FFT outdir)
#psd of fft (tsv and plot); plot excludes DC since data was mean-centred (saved in FFT outdir)
rule curvature_analyses_spaceUnfold:
    input:
        curvature_Unfolded_midthickness_nii = rules.curvatureT1w_to_curvatureUnfold_nii.output,
        label_subfields_spaceUnfold_nii = expand(rules.subfieldsT1w_to_spaceUnfold_nii.output, surface='midthickness', allow_missing=True),
    params:
        subfield_names = subfields
    output:
        curvature_subfields_analyses_Unfolded_to_smoothed_dir = directory('work/hippunfold_morphometry/sub-{subject}/curvature_analyses/{hemi}_surf-midthickness/curvature_mPD_at_eachAP/'),
        curvature_subfields_analyses_Unfolded_to_smoothed_FFTandPSD_dir = directory('work/hippunfold_morphometry/sub-{subject}/curvature_analyses/{hemi}_surf-midthickness/curvature_mPD_at_eachAP_FFTandPSD/'),
    group: 'subjLR'
    log: 'logs/hippunfold_morphometry/sub-{subject}/sub-{subject}_hemi-{hemi}_curvature_analyses_spaceUnfold.log'
    script:
        '../scripts/curvature_morphometry_midthickness_analysis.py'

#Group the tsv outputs of curvature_analyses_spaceUnfold into combined tsvs to permit future analyses on the group data
#    The output tsvs are created for the following subgroups: ipsilateral_TLE, contralateral_TLE, bilateral_TLE, allAffected_TLE (ipsilateral + bilateral), all_TLE (ipsilateral, contralateral, bilateral), controls. 
#      Headers have 2 columns: participant_id and <name_of_metric> (which is a list with length 87 of the  <name_of_metric> data for that subject)
#        <name_of_metric> can be one of the following values 'mPDcurvature_mAPcentered', 'FFT_mPDcurvature_mAPcentered', or 'PSDofFFT_mPDcurvature_mAPcentered'
#Currently only operating on CA1 and BodyOnly (see params)
rule group_curvature_analyses_spaceUnfold:
    input:
        curvature_subfields_analyses_Unfolded_to_smoothed_dir = expand(rules.curvature_analyses_spaceUnfold.output[0], subject=subjects, hemi=hemis, allow_missing=True),
        curvature_subfields_analyses_Unfolded_to_smoothed_FFTandPSD_dir = expand(rules.curvature_analyses_spaceUnfold.output[1], subject=subjects, hemi=hemis, allow_missing=True),
    params:
        demographics_tsv = config['demographics_tsv'],
        #Actual inputs
        L_CA1_tsv_subfield_curvature_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/curvature_analyses/L_surf-midthickness/curvature_mPD_at_eachAP/sub-{subject}_hemi-L_space-unfolded_subfield-CA1_BodyOnly_mPDcurvature_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
        L_CA1_tsv_FFT_subfield_curvature_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/curvature_analyses/L_surf-midthickness/curvature_mPD_at_eachAP_FFTandPSD/sub-{subject}_hemi-L_space-unfolded_subfield-CA1_BodyOnly_FFT_mPDcurvature_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
        L_CA1_tsv_FFT_PSD_subfield_curvature_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/curvature_analyses/L_surf-midthickness/curvature_mPD_at_eachAP_FFTandPSD/sub-{subject}_hemi-L_space-unfolded_subfield-CA1_BodyOnly_PSDofFFT_mPDcurvature_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
        R_CA1_tsv_subfield_curvature_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/curvature_analyses/R_surf-midthickness/curvature_mPD_at_eachAP/sub-{subject}_hemi-R_space-unfolded_subfield-CA1_BodyOnly_mPDcurvature_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
        R_CA1_tsv_FFT_subfield_curvature_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/curvature_analyses/R_surf-midthickness/curvature_mPD_at_eachAP_FFTandPSD/sub-{subject}_hemi-R_space-unfolded_subfield-CA1_BodyOnly_FFT_mPDcurvature_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
        R_CA1_tsv_FFT_PSD_subfield_curvature_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/curvature_analyses/R_surf-midthickness/curvature_mPD_at_eachAP_FFTandPSD/sub-{subject}_hemi-R_space-unfolded_subfield-CA1_BodyOnly_PSDofFFT_mPDcurvature_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
    output:
        group_CA1_tsv_subfield_curvature_spaceUnfold_mPD_atEachAP_dir = directory('results/hippunfold_morphometry/group_curvature_analyses/surface-midthickness/'),
    threads: 8
    group: 'groupSurfMid'
    log:
        'logs/hippunfold_morphometry/surf-midthickness_group_curvature_analyses.log'
    script:
        '../scripts/group_curvature_morphometry_surfs_analysis.py'

#Repeat for gyrification
rule gyrificationT1w_to_gyrificationUnfold_nii:
    input:
        gyrification_T1w = config['shape_gyrification_midthickness_spaceT1w'],
        gii_Unfolded_midthickness = expand(config['surf_spaceUnfold'], surface='midthickness', allow_missing=True),
        nii_Unfolded = config['refvol_spaceUnfold']
    params:
        gii_Unfolded_inner = expand(config['surf_spaceUnfold'], surface='inner', allow_missing=True),
        gii_Unfolded_outer = expand(config['surf_spaceUnfold'], surface='outer', allow_missing=True),
    output:
        gyrification_Unfolded_midthickness_nii = 'work/hippunfold_morphometry/sub-{subject}/anat/sub-{subject}_hemi-{hemi}_space-unfolded_den-0p5mm_label-hipp_midthickness_gyrification.nii.gz'
    group: 'subjLR'
    log: 'logs/hippunfold_morphometry/sub-{subject}/sub-{subject}_hemi-{hemi}_shape_gyrification_midthickness_T1w_xfm_spaceUnfold_nii.log'
    shell:
        'wb_command -metric-to-volume-mapping {input[0]} {input[1]} {input[2]} {output} -ribbon-constrained {params[0]} {params[1]} &> {log}'

rule gyrification_analyses_spaceUnfold:
    input:
        gyrification_Unfolded_midthickness_nii = rules.gyrificationT1w_to_gyrificationUnfold_nii.output,
        label_subfields_spaceUnfold_nii = expand(rules.subfieldsT1w_to_spaceUnfold_nii.output, surface='midthickness', allow_missing=True),
    params:
        subfield_names = subfields
    output:
        gyrification_subfields_analyses_Unfolded_to_smoothed_dir = directory('work/hippunfold_morphometry/sub-{subject}/gyrification_analyses/{hemi}_surf-midthickness/gyrification_mPD_at_eachAP/'),
        gyrification_subfields_analyses_Unfolded_to_smoothed_FFTandPSD_dir = directory('work/hippunfold_morphometry/sub-{subject}/gyrification_analyses/{hemi}_surf-midthickness/gyrification_mPD_at_eachAP_FFTandPSD/'),
    group: 'subjLR'
    log: 'logs/hippunfold_morphometry/sub-{subject}/sub-{subject}_hemi-{hemi}_gyrification_analyses_spaceUnfold.log'
    script:
        '../scripts/gyrification_morphometry_midthickness_analysis.py'

rule group_gyrification_analyses_spaceUnfold:
    input:
        gyrification_subfields_analyses_Unfolded_to_smoothed_dir = expand(rules.gyrification_analyses_spaceUnfold.output[0], subject=subjects, hemi=hemis, allow_missing=True),
        gyrification_subfields_analyses_Unfolded_to_smoothed_FFTandPSD_dir = expand(rules.gyrification_analyses_spaceUnfold.output[1], subject=subjects, hemi=hemis, allow_missing=True),
    params:
        demographics_tsv = config['demographics_tsv'],
        L_CA1_tsv_subfield_gyrification_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/gyrification_analyses/L_surf-midthickness/gyrification_mPD_at_eachAP/sub-{subject}_hemi-L_space-unfolded_subfield-CA1_BodyOnly_mPDgyrification_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
        L_CA1_tsv_FFT_subfield_gyrification_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/gyrification_analyses/L_surf-midthickness/gyrification_mPD_at_eachAP_FFTandPSD/sub-{subject}_hemi-L_space-unfolded_subfield-CA1_BodyOnly_FFT_mPDgyrification_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
        L_CA1_tsv_FFT_PSD_subfield_gyrification_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/gyrification_analyses/L_surf-midthickness/gyrification_mPD_at_eachAP_FFTandPSD/sub-{subject}_hemi-L_space-unfolded_subfield-CA1_BodyOnly_PSDofFFT_mPDgyrification_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
        R_CA1_tsv_subfield_gyrification_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/gyrification_analyses/R_surf-midthickness/gyrification_mPD_at_eachAP/sub-{subject}_hemi-R_space-unfolded_subfield-CA1_BodyOnly_mPDgyrification_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
        R_CA1_tsv_FFT_subfield_gyrification_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/gyrification_analyses/R_surf-midthickness/gyrification_mPD_at_eachAP_FFTandPSD/sub-{subject}_hemi-R_space-unfolded_subfield-CA1_BodyOnly_FFT_mPDgyrification_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
        R_CA1_tsv_FFT_PSD_subfield_gyrification_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/gyrification_analyses/R_surf-midthickness/gyrification_mPD_at_eachAP_FFTandPSD/sub-{subject}_hemi-R_space-unfolded_subfield-CA1_BodyOnly_PSDofFFT_mPDgyrification_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
    output:
        group_CA1_tsv_subfield_gyrification_spaceUnfold_mPD_atEachAP_dir = directory('results/hippunfold_morphometry/group_gyrification_analyses/surface-midthickness/'),
    threads: 8
    group: 'groupSurfMid'
    log:
        'logs/hippunfold_morphometry/surf-midthickness_group_gyrification_analyses.log'
    script:
        '../scripts/group_gyrification_morphometry_surfs_analysis.py'
