#Thickness
#Map thickness to unfolded VOLUME space, using unfolded midsurf and constrained by unfolded inner/outer surfs
rule thicknessT1w_to_thicknessUnfold_nii:
    input:
        thickness_T1w = config['shape_thickness_midthickness_spaceT1w'], #Metric to map to vol
        gii_Unfolded_midthickness = expand(config['surf_spaceUnfold'], surface='midthickness', allow_missing=True), #Unfolded midsurf
        nii_Unfolded = config['refvol_spaceUnfold'] #Ref Unfolded vol to map metric to
    params:
        gii_Unfolded_inner = expand(config['surf_spaceUnfold'], surface='inner', allow_missing=True), #For -ribbon-constrained
        gii_Unfolded_outer = expand(config['surf_spaceUnfold'], surface='outer', allow_missing=True), #For -ribbon-constrained
    output:
        thickness_Unfolded_midthickness_nii = 'work/hippunfold_morphometry/sub-{subject}/anat/sub-{subject}_hemi-{hemi}_space-unfolded_den-0p5mm_label-hipp_midthickness_thickness.nii.gz'
    group: 'subjLR'
    log: 'logs/hippunfold_morphometry/sub-{subject}/sub-{subject}_hemi-{hemi}_shape_thickness_midthickness_T1w_xfm_spaceUnfold_nii.log'
    shell:
        'wb_command -metric-to-volume-mapping {input[0]} {input[1]} {input[2]} {output} -ribbon-constrained {params[0]} {params[1]} &> {log}'


#Thickness analyses, outputs are as follows and all exist for each subfield and split into wholeHc(5% clipped from head and tail), HeadOnly, and BodyOnly (individual surfaces and hemis are processed in separate calls).
#Each AP's average (across PD) thickness values (mean-centred) in each individual subfield (tsv and plot)
#fft of this (tsv and plot); plot excludes DC since data was mean-centred (saved in FFT outdir)
#psd of fft (tsv and plot); plot excludes DC since data was mean-centred (saved in FFT outdir)
rule thickness_analyses_spaceUnfold:
    input:
        thickness_Unfolded_midthickness_nii = rules.thicknessT1w_to_thicknessUnfold_nii.output,
        label_subfields_spaceUnfold_nii = expand(rules.subfieldsT1w_to_spaceUnfold_nii.output, surface='midthickness', allow_missing=True),
    params:
        subfield_names = subfields
    output:
        thickness_subfields_analyses_Unfolded_to_smoothed_dir = directory('work/hippunfold_morphometry/sub-{subject}/thickness_analyses/{hemi}_surf-midthickness/thickness_mPD_at_eachAP/'),
        thickness_subfields_analyses_Unfolded_to_smoothed_FFTandPSD_dir = directory('work/hippunfold_morphometry/sub-{subject}/thickness_analyses/{hemi}_surf-midthickness/thickness_mPD_at_eachAP_FFTandPSD/'),
    group: 'subjLR'
    log: 'logs/hippunfold_morphometry/sub-{subject}/sub-{subject}_hemi-{hemi}_thickness_analyses_spaceUnfold.log'
    script:
        '../scripts/thickness_morphometry_midthickness_analysis.py'

#Group the tsv outputs of thickness_analyses_spaceUnfold into combined tsvs to permit future analyses on the group data
#    The output tsvs are created for the following subgroups: ipsilateral_TLE, contralateral_TLE, bilateral_TLE, allAffected_TLE (ipsilateral + bilateral), all_TLE (ipsilateral, contralateral, bilateral), controls. 
#      Headers have 2 columns: participant_id and <name_of_metric> (which is a list with length 87 of the  <name_of_metric> data for that subject)
#        <name_of_metric> can be one of the following values 'mPDthickness_mAPcentered', 'FFT_mPDthickness_mAPcentered', or 'PSDofFFT_mPDthickness_mAPcentered'
#Currently only operating on CA1 and BodyOnly (see params)
rule group_thickness_analyses_spaceUnfold:
    input:
        thickness_subfields_analyses_Unfolded_to_smoothed_dir = expand(rules.thickness_analyses_spaceUnfold.output[0], subject=subjects, hemi=hemis, allow_missing=True),
        thickness_subfields_analyses_Unfolded_to_smoothed_FFTandPSD_dir = expand(rules.thickness_analyses_spaceUnfold.output[1], subject=subjects, hemi=hemis, allow_missing=True),
    params:
        demographics_tsv = config['demographics_tsv'],
        #Actual inputs
        L_CA1_tsv_subfield_thickness_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/thickness_analyses/L_surf-midthickness/thickness_mPD_at_eachAP/sub-{subject}_hemi-L_space-unfolded_subfield-CA1_BodyOnly_mPDthickness_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
        L_CA1_tsv_FFT_subfield_thickness_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/thickness_analyses/L_surf-midthickness/thickness_mPD_at_eachAP_FFTandPSD/sub-{subject}_hemi-L_space-unfolded_subfield-CA1_BodyOnly_FFT_mPDthickness_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
        L_CA1_tsv_FFT_PSD_subfield_thickness_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/thickness_analyses/L_surf-midthickness/thickness_mPD_at_eachAP_FFTandPSD/sub-{subject}_hemi-L_space-unfolded_subfield-CA1_BodyOnly_PSDofFFT_mPDthickness_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
        R_CA1_tsv_subfield_thickness_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/thickness_analyses/R_surf-midthickness/thickness_mPD_at_eachAP/sub-{subject}_hemi-R_space-unfolded_subfield-CA1_BodyOnly_mPDthickness_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
        R_CA1_tsv_FFT_subfield_thickness_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/thickness_analyses/R_surf-midthickness/thickness_mPD_at_eachAP_FFTandPSD/sub-{subject}_hemi-R_space-unfolded_subfield-CA1_BodyOnly_FFT_mPDthickness_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
        R_CA1_tsv_FFT_PSD_subfield_thickness_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/thickness_analyses/R_surf-midthickness/thickness_mPD_at_eachAP_FFTandPSD/sub-{subject}_hemi-R_space-unfolded_subfield-CA1_BodyOnly_PSDofFFT_mPDthickness_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
    output:
        group_CA1_tsv_subfield_thickness_spaceUnfold_mPD_atEachAP_dir = directory('results/hippunfold_morphometry/group_thickness_analyses/surface-midthickness/'),
    threads: 8
    group: 'groupSurfMid'
    log:
        'logs/hippunfold_morphometry/surf-midthickness_group_thickness_analyses.log'
    script:
        '../scripts/group_thickness_morphometry_surfs_analysis.py'

#Repeat for surfarea
rule surfareaT1w_to_surfareaUnfold_nii:
    input:
        surfarea_T1w = config['shape_surfarea_midthickness_spaceT1w'],
        gii_Unfolded_midthickness = expand(config['surf_spaceUnfold'], surface='midthickness', allow_missing=True),
        nii_Unfolded = config['refvol_spaceUnfold']
    params:
        gii_Unfolded_inner = expand(config['surf_spaceUnfold'], surface='inner', allow_missing=True),
        gii_Unfolded_outer = expand(config['surf_spaceUnfold'], surface='outer', allow_missing=True),
    output:
        surfarea_Unfolded_midthickness_nii = 'work/hippunfold_morphometry/sub-{subject}/anat/sub-{subject}_hemi-{hemi}_space-unfolded_den-0p5mm_label-hipp_midthickness_surfarea.nii.gz'
    group: 'subjLR'
    log: 'logs/hippunfold_morphometry/sub-{subject}/sub-{subject}_hemi-{hemi}_shape_surfarea_midthickness_T1w_xfm_spaceUnfold_nii.log'
    shell:
        'wb_command -metric-to-volume-mapping {input[0]} {input[1]} {input[2]} {output} -ribbon-constrained {params[0]} {params[1]} &> {log}'

rule surfarea_analyses_spaceUnfold:
    input:
        surfarea_Unfolded_midthickness_nii = rules.surfareaT1w_to_surfareaUnfold_nii.output,
        label_subfields_spaceUnfold_nii = expand(rules.subfieldsT1w_to_spaceUnfold_nii.output, surface='midthickness', allow_missing=True),
    params:
        subfield_names = subfields
    output:
        surfarea_subfields_analyses_Unfolded_to_smoothed_dir = directory('work/hippunfold_morphometry/sub-{subject}/surfarea_analyses/{hemi}_surf-midthickness/surfarea_mPD_at_eachAP/'),
        surfarea_subfields_analyses_Unfolded_to_smoothed_FFTandPSD_dir = directory('work/hippunfold_morphometry/sub-{subject}/surfarea_analyses/{hemi}_surf-midthickness/surfarea_mPD_at_eachAP_FFTandPSD/'),
    group: 'subjLR'
    log: 'logs/hippunfold_morphometry/sub-{subject}/sub-{subject}_hemi-{hemi}_surfarea_analyses_spaceUnfold.log'
    script:
        '../scripts/surfarea_morphometry_midthickness_analysis.py'

rule group_surfarea_analyses_spaceUnfold:
    input:
        surfarea_subfields_analyses_Unfolded_to_smoothed_dir = expand(rules.surfarea_analyses_spaceUnfold.output[0], subject=subjects, hemi=hemis, allow_missing=True),
        surfarea_subfields_analyses_Unfolded_to_smoothed_FFTandPSD_dir = expand(rules.surfarea_analyses_spaceUnfold.output[1], subject=subjects, hemi=hemis, allow_missing=True),
    params:
        demographics_tsv = config['demographics_tsv'],
        L_CA1_tsv_subfield_surfarea_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/surfarea_analyses/L_surf-midthickness/surfarea_mPD_at_eachAP/sub-{subject}_hemi-L_space-unfolded_subfield-CA1_BodyOnly_mPDsurfarea_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
        L_CA1_tsv_FFT_subfield_surfarea_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/surfarea_analyses/L_surf-midthickness/surfarea_mPD_at_eachAP_FFTandPSD/sub-{subject}_hemi-L_space-unfolded_subfield-CA1_BodyOnly_FFT_mPDsurfarea_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
        L_CA1_tsv_FFT_PSD_subfield_surfarea_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/surfarea_analyses/L_surf-midthickness/surfarea_mPD_at_eachAP_FFTandPSD/sub-{subject}_hemi-L_space-unfolded_subfield-CA1_BodyOnly_PSDofFFT_mPDsurfarea_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
        R_CA1_tsv_subfield_surfarea_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/surfarea_analyses/R_surf-midthickness/surfarea_mPD_at_eachAP/sub-{subject}_hemi-R_space-unfolded_subfield-CA1_BodyOnly_mPDsurfarea_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
        R_CA1_tsv_FFT_subfield_surfarea_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/surfarea_analyses/R_surf-midthickness/surfarea_mPD_at_eachAP_FFTandPSD/sub-{subject}_hemi-R_space-unfolded_subfield-CA1_BodyOnly_FFT_mPDsurfarea_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
        R_CA1_tsv_FFT_PSD_subfield_surfarea_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/surfarea_analyses/R_surf-midthickness/surfarea_mPD_at_eachAP_FFTandPSD/sub-{subject}_hemi-R_space-unfolded_subfield-CA1_BodyOnly_PSDofFFT_mPDsurfarea_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
    output:
        group_CA1_tsv_subfield_surfarea_spaceUnfold_mPD_atEachAP_dir = directory('results/hippunfold_morphometry/group_surfarea_analyses/surface-midthickness/'),
    threads: 8
    group: 'groupSurfMid'
    log:
        'logs/hippunfold_morphometry/surf-midthickness_group_surfarea_analyses.log'
    script:
        '../scripts/group_surfarea_morphometry_surfs_analysis.py'
