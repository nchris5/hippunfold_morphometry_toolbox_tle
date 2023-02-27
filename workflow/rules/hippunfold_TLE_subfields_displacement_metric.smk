#Displacement

#Smooth T1w surfaces (inner, midthickness, outer) ... tested various smoothing strengths(0.5 0.7 0.9) and iterations = 10 20 30 50 100. Seems best with 0.9 and 50+ iterations so this is what is set in the config
#  To test other values change surf_smoothing_strength/iter in the config and use the following while editing in vim 
#    :%s/str<current_smooth_strength>/str<new_smooth strength>/g or :%s/iter<current_iter>/iter<new_iter>/g
rule surf_smooth:
    input:
        gii_T1w = config['surf_spaceT1w'],
    params:
        surf_smooth_s = config['surf_smoothing_strength'],
        surf_smooth_i = config['surf_smoothing_iter'],
    output:
        gii_T1w_smoothed = 'work/hippunfold_morphometry/sub-{subject}/surf/sub-{subject}_hemi-{hemi}_surface-{surface}_smooth/sub-{subject}_hemi-{hemi}_space-T1w_den-0p5mm_label-hipp_{surface}_smoothed_str0.9_iter50.surf.gii',
    log: 'logs/hippunfold_morphometry/sub-{subject}/sub-{subject}_hemi-{hemi}_surf-{surface}_smooth.log'
    threads: 8
    group: 'surfsLR'
    shell:
        'wb_command -surface-smoothing {input} {params[0]} {params[1]} {output} &> {log}'

#Lightly smooth for downstream displacement calculation vs. fully smoothed surface, this removes high-frequency variations, while still maintaining low-frequency digitations
rule surf_smoothSoft:
    input:
        gii_T1w = config['surf_spaceT1w'],
    params:
        surf_smoothSoft_s = config['surf_smoothing_soft_strength'],
        surf_smoothSoft_i = config['surf_smoothing_soft_iter'],
    output:
        gii_T1w_smoothedSoft = 'work/hippunfold_morphometry/sub-{subject}/surf/sub-{subject}_hemi-{hemi}_surface-{surface}_smooth/sub-{subject}_hemi-{hemi}_space-T1w_den-0p5mm_label-hipp_{surface}_smoothed_str0.5_iter10.surf.gii',
    log: 'logs/hippunfold_morphometry/sub-{subject}/sub-{subject}_hemi-{hemi}_surf-{surface}_smooth_soft.log'
    threads: 8
    group: 'surfsLR'
    shell:
        'wb_command -surface-smoothing {input} {params[0]} {params[1]} {output} &> {log}'


#Find displacement from surfs to smoothed_surfs (inner midthickness outer) (consider using -vectors to save displacement vectors)
#This is done for both spaceT1w displacementTo spaceT1w_fullySmoothed ... and
#                      spaceT1w_softSmooted displacementTo spaceT1w_fullySmoothed
rule surf_to_smoothSurf_displacement:
    input:
        gii_T1w = config['surf_spaceT1w'], #Input unsmoothed
        gii_T1w_smoothedSoft = rules.surf_smoothSoft.output, #Input smoothedSoft
        gii_T1w_smoothed = rules.surf_smooth.output, #Reference smoothedHard
    output:
        displacement_T1w_to_smoothed = 'work/hippunfold_morphometry/sub-{subject}/surf/sub-{subject}_hemi-{hemi}_space-T1w_den-0p5mm_label-hipp_{surface}_displacement_to_smoothed_str0.9_iter50.shape.gii',
        displacement_T1w_smoothedSoft_to_smoothed = 'work/hippunfold_morphometry/sub-{subject}/surf/sub-{subject}_hemi-{hemi}_space-T1w_smoothedSoft_den-0p5mm_label-hipp_{surface}_displacement_to_smoothed_str0.9_iter50.shape.gii',
    log: 'logs/hippunfold_morphometry/sub-{subject}/sub-{subject}_hemi-{hemi}_surf-{surface}_surf_to_smoothSurf_displacement.log'
    threads: 8
    group: 'surfsLR'
    shell:
        'wb_command -surface-to-surface-3d-distance {input[0]} {input[2]} {output[0]} && '
        'wb_command -surface-to-surface-3d-distance {input[1]} {input[2]} {output[1]} &> {log}'


###Map displacement to unfolded VOLUME space, using unfolded surfs (inner midthickness outer) and constrained by unfolded inner/outer surfs
#This is done for both spaceT1w displacementTo spaceT1w_fullySmoothed_nii ... and
#                      spaceT1w_softSmooted displacementTo spaceT1w_fullySmoothed_nii
rule displacementT1w_to_displacementUnfold_nii:
    input:
        displacement_T1w_to_smoothed = rules.surf_to_smoothSurf_displacement.output[0], #Metric to map to vol
        displacement_T1w_smoothedSoft_to_smoothed = rules.surf_to_smoothSurf_displacement.output[1],
        gii_Unfolded = config['surf_spaceUnfold'], #Unfolded surfs (inner mid outer)
        nii_Unfolded = config['refvol_spaceUnfold'] #Ref Unfolded vol to map metric to
    params:
        gii_Unfolded_inner = expand(config['surf_spaceUnfold'], surface='inner', allow_missing=True), #For -ribbon-constrained
        gii_Unfolded_outer = expand(config['surf_spaceUnfold'], surface='outer', allow_missing=True), #For -ribbon-constrained
    output:
        displacement_Unfolded_to_smoothed_nii = 'work/hippunfold_morphometry/sub-{subject}/anat/sub-{subject}_hemi-{hemi}_space-unfolded_den-0p5mm_label-hipp_{surface}_displacement_to_smoothed_str0.9_iter50.nii.gz',
        displacement_Unfolded_smoothedSoft_to_smoothed_nii = 'work/hippunfold_morphometry/sub-{subject}/anat/sub-{subject}_hemi-{hemi}_space-unfolded_den-0p5mm_label-hipp_{surface}_soft_smoothed_displacement_to_smoothed_str0.9_iter50.nii.gz',
    group: 'surfsLR'
    log: 'logs/hippunfold_morphometry/sub-{subject}/sub-{subject}_hemi-{hemi}_{surface}_shape_displacement_T1w_xfm_spaceUnfold_nii.log'
    threads: 8
    shell:
        'wb_command -metric-to-volume-mapping {input[0]} {input[2]} {input[3]} {output[0]} -ribbon-constrained {params[0]} {params[1]} && '
        'wb_command -metric-to-volume-mapping {input[1]} {input[2]} {input[3]} {output[1]} -ribbon-constrained {params[0]} {params[1]} &> {log}'


#Same for subfield labels (this was similarly done in previous rule)
#Only done for unsmoothed surface
rule subfieldsT1w_to_spaceUnfold_nii:
    input:
        label_subfieldsT1w = config['label_subfields_spaceT1w'], #Metric to map to vol
        gii_Unfolded = config['surf_spaceUnfold'], #Unfolded surfs (inner mid outer)
        nii_Unfolded = config['refvol_spaceUnfold'] #Ref vol to map metric to
    params:
        gii_Unfolded_inner = expand(config['surf_spaceUnfold'], surface='inner', allow_missing=True), #For -ribbon-constrained
        gii_Unfolded_outer = expand(config['surf_spaceUnfold'], surface='outer', allow_missing=True), #For -ribbon-constrained
    output:
        labelSubfields_Unfolded_nii = 'work/hippunfold_morphometry/sub-{subject}/anat/sub-{subject}_hemi-{hemi}_space-unfolded_den-0p5mm_label-hipp_{surface}_atlas-bigbrain_subfields.nii.gz',
    group: 'surfsLR'
    log: 'logs/hippunfold_morphometry/sub-{subject}/sub-{subject}_hemi-{hemi}_{surface}_subfieldlabels_surfs_T1w_xfm_spaceUnfold_nii.log'
    threads: 8
    shell:
        'wb_command -label-to-volume-mapping {input[0]} {input[1]} {input[2]} {output} -ribbon-constrained {params[0]} {params[1]} &> {log}'


#Displacement analyses, outputs are as follows and all exist for each subfield and split into wholeHc(5% clipped from head and tail), HeadOnly, and BodyOnly (individual surfaces and hemis are processed in separate calls).
#Each AP's average (across PD) displacement values (mean-centred) in each individual subfield (tsv and plot)
#fft of this (tsv and plot); plot excludes DC since data was mean-centred (saved in FFT outdir)
#psd of fft (tsv and plot); plot excludes DC since data was mean-centred (saved in FFT outdir)
rule displacement_analyses_spaceUnfold:
    input:
        displacement_Unfolded_to_smoothed_nii = rules.displacementT1w_to_displacementUnfold_nii.output[0],
        labelSubfields_Unfolded_nii = rules.subfieldsT1w_to_spaceUnfold_nii.output[0],
    params:
        subfield_names = subfields
    output:
        displacement_subfields_analyses_Unfolded_to_smoothed_dir = directory('work/hippunfold_morphometry/sub-{subject}/displacement_analysis/{hemi}_surf-{surface}/displacement_mPD_at_eachAP/'),
        displacement_subfields_analyses_Unfolded_to_smoothed_FFTandPSD_dir = directory('work/hippunfold_morphometry/sub-{subject}/displacement_analysis/{hemi}_surf-{surface}/displacement_mPD_at_eachAP_FFTandPSD/'),
    group: 'surfsLR'
    log: 'logs/hippunfold_morphometry/sub-{subject}/sub-{subject}_hemi-{hemi}_{surface}_displacement_analyses_spaceUnfold.log'
    threads: 8
    script:
        '../scripts/displacement_morphometry_surfs_analysis.py'

#Group the tsv outputs of displacement_analyses_spaceUnfold into combined tsvs to permit future analyses on the group data
#    The output tsvs are created for the following subgroups: ipsilateral_TLE, contralateral_TLE, bilateral_TLE, allAffected_TLE (ipsilateral + bilateral), all_TLE (ipsilateral, contralateral, bilateral), controls. 
#      Headers have 2 columns: participant_id and <name_of_metric> (which is a list with length 87 of the  <name_of_metric> data for that subject)
#        <name_of_metric> can be one of the following values 'mPDdisplacement_mAPcentered', 'FFT_mPDdisplacement_mAPcentered', or 'PSDofFFT_mPDdisplacement_mAPcentered'
#Currently only operating on CA1 and BodyOnly (see params)
rule group_displacement_analyses_spaceUnfold:
    input:
        displacement_subfields_analyses_Unfolded_to_smoothed_dir = expand(rules.displacement_analyses_spaceUnfold.output[0], subject=subjects, hemi=hemis, allow_missing=True),
        displacement_subfields_analyses_Unfolded_to_smoothed_FFTandPSD_dir = expand(rules.displacement_analyses_spaceUnfold.output[1], subject=subjects, hemi=hemis, allow_missing=True),
    params:
        demographics_tsv = config['demographics_tsv'],
        #Actual inputs (indexing inputs in python starts at 1)
        L_CA1_tsv_subfield_displacement_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/displacement_analysis/L_surf-{surface}/displacement_mPD_at_eachAP/sub-{subject}_hemi-L_space-unfolded_subfield-CA1_BodyOnly_mPDdisplacement_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
        L_CA1_tsv_FFT_subfield_displacement_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/displacement_analysis/L_surf-{surface}/displacement_mPD_at_eachAP_FFTandPSD/sub-{subject}_hemi-L_space-unfolded_subfield-CA1_BodyOnly_FFT_mPDdisplacement_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
        L_CA1_tsv_FFT_PSD_subfield_displacement_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/displacement_analysis/L_surf-{surface}/displacement_mPD_at_eachAP_FFTandPSD/sub-{subject}_hemi-L_space-unfolded_subfield-CA1_BodyOnly_PSDofFFT_mPDdisplacement_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
        R_CA1_tsv_subfield_displacement_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/displacement_analysis/R_surf-{surface}/displacement_mPD_at_eachAP/sub-{subject}_hemi-R_space-unfolded_subfield-CA1_BodyOnly_mPDdisplacement_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
        R_CA1_tsv_FFT_subfield_displacement_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/displacement_analysis/R_surf-{surface}/displacement_mPD_at_eachAP_FFTandPSD/sub-{subject}_hemi-R_space-unfolded_subfield-CA1_BodyOnly_FFT_mPDdisplacement_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
        R_CA1_tsv_FFT_PSD_subfield_displacement_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/displacement_analysis/R_surf-{surface}/displacement_mPD_at_eachAP_FFTandPSD/sub-{subject}_hemi-R_space-unfolded_subfield-CA1_BodyOnly_PSDofFFT_mPDdisplacement_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
    output:
        group_CA1_tsv_subfield_displacement_spaceUnfold_mPD_atEachAP_dir = directory('results/hippunfold_morphometry/group_displacement_analysis/surface-{surface}/'),
    threads: 8
    group: 'groupSurfs'
    log:
        'logs/hippunfold_morphometry/surf-{surface}_group_displacement_analyses.log'
    script:
        '../scripts/group_displacement_morphometry_surfs_analysis.py'

#Repeat displacement analyses and grouping for smoothedSoft_spaceUnfold to smoothed_spaceUnfold
rule displacement_analyses_smoothedSoft_spaceUnfold:
    input:
        displacement_T1w_smoothedSoft_to_smoothed = rules.displacementT1w_to_displacementUnfold_nii.output[1],
        label_subfields_spaceUnfold_nii = rules.subfieldsT1w_to_spaceUnfold_nii.output[0],
    params:
        subfield_names = subfields
    output:
        displacement_subfields_analyses_softSmoothed_Unfolded_to_smoothed_dir = directory('work/hippunfold_morphometry/sub-{subject}/displacement_analysis_smoothedSoft/{hemi}_surf-{surface}/displacement_mPD_at_eachAP/'),
        displacement_subfields_analyses_softSmoothed_Unfolded_to_smoothed_FFTandPSD_dir = directory('work/hippunfold_morphometry/sub-{subject}/displacement_analysis_smoothedSoft/{hemi}_surf-{surface}/displacement_mPD_at_eachAP_FFTandPSD/'),
    group: 'surfsLR'
    log: 'logs/hippunfold_morphometry/sub-{subject}/sub-{subject}_hemi-{hemi}_{surface}_displacement_analyses_smoothed_soft_spaceUnfold.log'
    threads: 8
    script:
        '../scripts/displacement_morphometry_surfs_analysis.py'

rule group_displacement_analyses_smoothedSoft_spaceUnfold:
    input:
        displacement_subfields_analyses_softSmoothed_Unfolded_to_smoothed_dir = expand(rules.displacement_analyses_smoothedSoft_spaceUnfold.output[0], subject=subjects, hemi=hemis, allow_missing=True),
        displacement_subfields_analyses_softSmoothed_Unfolded_to_smoothed_FFTandPSD_dir = expand(rules.displacement_analyses_smoothedSoft_spaceUnfold.output[1], subject=subjects, hemi=hemis, allow_missing=True),
    params:
        demographics_tsv = config['demographics_tsv'],
        #Actual inputs
        L_CA1_tsv_subfield_displacement_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/displacement_analysis_smoothedSoft/L_surf-{surface}/displacement_mPD_at_eachAP/sub-{subject}_hemi-L_space-unfolded_subfield-CA1_BodyOnly_mPDdisplacement_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
        L_CA1_tsv_FFT_subfield_displacement_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/displacement_analysis_smoothedSoft/L_surf-{surface}/displacement_mPD_at_eachAP_FFTandPSD/sub-{subject}_hemi-L_space-unfolded_subfield-CA1_BodyOnly_FFT_mPDdisplacement_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
        L_CA1_tsv_FFT_PSD_subfield_displacement_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/displacement_analysis_smoothedSoft/L_surf-{surface}/displacement_mPD_at_eachAP_FFTandPSD/sub-{subject}_hemi-L_space-unfolded_subfield-CA1_BodyOnly_PSDofFFT_mPDdisplacement_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
        R_CA1_tsv_subfield_displacement_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/displacement_analysis_smoothedSoft/R_surf-{surface}/displacement_mPD_at_eachAP/sub-{subject}_hemi-R_space-unfolded_subfield-CA1_BodyOnly_mPDdisplacement_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
        R_CA1_tsv_FFT_subfield_displacement_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/displacement_analysis_smoothedSoft/R_surf-{surface}/displacement_mPD_at_eachAP_FFTandPSD/sub-{subject}_hemi-R_space-unfolded_subfield-CA1_BodyOnly_FFT_mPDdisplacement_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
        R_CA1_tsv_FFT_PSD_subfield_displacement_spaceUnfold_mPD_atEachAP = expand('work/hippunfold_morphometry/sub-{subject}/displacement_analysis_smoothedSoft/R_surf-{surface}/displacement_mPD_at_eachAP_FFTandPSD/sub-{subject}_hemi-R_space-unfolded_subfield-CA1_BodyOnly_PSDofFFT_mPDdisplacement_mAPcentered_at_each_AP.tsv', subject=subjects, allow_missing=True),
    output:
        group_CA1_tsv_subfield_displacement_spaceUnfold_mPD_atEachAP_dir = directory('results/hippunfold_morphometry/group_displacement_analysis_smoothedSoft/surface-{surface}/'),
    group: 'groupSurfs'
    log:
        'logs/hippunfold_morphometry/surf-{surface}_group_displacement_analyses_smoothed_soft.log'
    threads: 8
    script:
        '../scripts/group_displacement_morphometry_surfs_analysis.py'



#-------------- IN DEV --------------#
#Descriptive group statistics (Mean Std, Sterror) of rule  group_displacement_analyses_spaceUnfold output displacement values, FFT, and PSD.
#This separates groups based on ipsilateral_TLE, contralateral_TLE, controls, bilateral_TLE, ipsilateral+bilateral_TLE (allAffected_TLE), and ipsilateral+contralateral+bilateral_TLE (all_TLE)
#The following outputs are returned:
#1)  Each groups individual displacement values, FFT, and PSD for each subject (tsvs, each shape is #subjects * metric placed in a list in the 2nd column, which has size 87)
#2)  Each groups individual mean/std/sterror of displacement values, FFT, and PSD across subjects (tsvs, each shape is length of metric (87) * 3 (cols = mean, std, sterror)
#3)  Line plots of mean/sterrors for each group (single graph plots all groups):
#rule stats_group_displacement_analyses_spaceUnfold_CA1:
#    input:
#        group_CA1_tsv_subfield_displacement_spaceUnfold_mPD_atEachAP_dir = rules.group_displacement_analyses_spaceUnfold.output.group_CA1_tsv_subfield_displacement_spaceUnfold_mPD_atEachAP_dir,
#    params:
#        #Actual inputs (displacement values, FFT of displacement, PSD of displacement)
#        group_LRcombined_CA1_tsv_subfield_displacement_spaceUnfold_mPD_atEachAP = 'work/hippunfold_morphometry/group_displacement_analysis/surface-{surface}/group_hemi-LR_label-hipp_{surface}_CA1_BodyOnly_mPDdisplacement_at_each_AP.tsv',
#        group_LRcombined_CA1_tsv_subfield_FFT_displacement_spaceUnfold_mPD_atEachAP = 'work/hippunfold_morphometry/group_displacement_analysis/surface-{surface}/group_hemi-LR_label-hipp_{surface}_CA1_BodyOnly_FFT_mPDdisplacement_at_each_AP.tsv',
#        group_LRcombined_CA1_tsv_subfield_FFT_PSD_displacement_spaceUnfold_mPD_atEachAP = 'work/hippunfold_morphometry/group_displacement_analysis/surface-{surface}/group_hemi-LR_label-hipp_{surface}_CA1_BodyOnly_FFT_PSD_mPDdisplacement_at_each_AP.tsv',
#    output:
#        stats_group_LRcombined_CA1_tsv_subfield_displacement_spaceUnfold_mPD_atEachAP_dir = directory('work/hippunfold_morphometry/group_displacement_analysis/surface-{surface}_stats/displacement_value_analysis/'),
#        stats_group_LRcombined_CA1_tsv_subfield_FFT_displacement_spaceUnfold_mPD_atEachAP_dir = directory('work/hippunfold_morphometry/group_displacement_analysis/surface-{surface}_stats/displacement_value_FFT_analysis/'),
#        stats_group_LRcombined_CA1_tsv_subfield_FFT_PSD_displacement_spaceUnfold_mPD_atEachAP_dir = directory('work/hippunfold_morphometry/group_displacement_analysis/surface-{surface}_stats/displacement_value_FFT_PSD_analysis/'),
#    group: 'groupSurfs'
#    log:
#        'logs/hippunfold_morphometry/stats_surf-{surface}_group_displacement_analyses.log'
#    threads: 8
#    script:
#        '../scripts/stats_group_displacement_morphometry_surfs_analysis.py'
#    ###For some reason PSD comes out with size 2 smaller instead of expected, need to figure out why this is
