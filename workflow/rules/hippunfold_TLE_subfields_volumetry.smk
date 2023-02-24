#Volumetry
rule volumetry_T1w_subfields:
    input:
        dseg_T1w_subfields = config['dseg_subfields_bigbrain_spaceT1w']
    params:
        voxel_shape = config['voxel_shape']
    output:
        volumetry_T1w_subfields = 'work/hippunfold_volumetry/sub-{subject}/sub-{subject}_hemi-{hemi}_space-T1w_desc-subfields_atlas-bigbrain_volumes.tsv'
    log: 'logs/hippunfold_volumetry/sub-{subject}_hemi-{hemi}_volumetry_T1w_subfields_dseg.log'
    group: 'volLR'
    script: '../scripts/volumetry_subfields.py'

#Group individual subjects hippocampal subfields into a single tsv file
rule group_volumetry_T1w_subfields:
    input:
        L_volumetry_T1w_subfields_tsv = expand(rules.volumetry_T1w_subfields.output, subject=subjects, hemi='L'),
        R_volumetry_T1w_subfields_tsv = expand(rules.volumetry_T1w_subfields.output, subject=subjects, hemi='R'),
    params:
        demographics_tsv = config['demographics_tsv']
    output:
        group_volumetry_T1w_subfields_tsv_dir = directory('results/hippunfold_volumetry/group_volumetry/')
    log: 'logs/hippunfold_volumetry/group_volumetry_T1w_subfields_dseg.log'
    group: 'groupVol'
    script: '../scripts/group_volumetry_subfields.py'

###-------------INDEV--------------
#rule tTest_group_hippunfold_subfields_T1w_dseg_volumetry:
#    input:
#        subfield_volumes_subj_group_tsv = 'work/subfield_volumes/group_space-T1w_desc-subfields_atlas-bigbrain_dseg_volumetry.tsv'
#    output:
#        tTest_stats_wholeHc_TLEips_TLEcon_tsv = 'work/subfield_volumes_statistics/tTest_stats_wholeHc_TLEips_TLEcon_space-T1w_desc-subfields_atlas-bigbrain_dseg_volumetry.tsv',
#        tTest_stats_wholeHc_TLEips_Control_tsv = 'work/subfield_volumes_statistics/tTest_stats_wholeHc_TLEips_Control_space-T1w_desc-subfields_atlas-bigbrain_dseg_volumetry.tsv',
#        tTest_stats_subfields_TLEips_TLEcon_tsv = 'work/subfield_volumes_statistics/tTest_stats_subfieldsHc_TLEips_TLEcon_space-T1w_desc-subfields_atlas-bigbrain_dseg_volumetry.tsv',
#        tTest_stats_subfields_TLEips_Control_tsv = 'work/subfield_volumes_statistics/tTest_stats_subfieldsHc_TLEips_Control_space-T1w_desc-subfields_atlas-bigbrain_dseg_volumetry.tsv',
#        tTest_stats_wholeHc_TLEips_TLEcon_fig = 'work/subfield_volumes_statistics/tTest_stats_wholeHc_TLEips_TLEcon_space-T1w_desc-subfields_atlas-bigbrain_dseg_volumetry.png',
#        tTest_stats_wholeHc_TLEips_Control_fig = 'work/subfield_volumes_statistics/tTest_stats_wholeHc_TLEips_Control_space-T1w_desc-subfields_atlas-bigbrain_dseg_volumetry.png',
#        tTest_stats_subfields_TLEips_TLEcon_fig = 'work/subfield_volumes_statistics/tTest_stats_subfieldsHc_TLEips_TLEcon_space-T1w_desc-subfields_atlas-bigbrain_dseg_volumetry.png',
#        tTest_stats_subfields_TLEips_Control_fig = 'work/subfield_volumes_statistics/tTest_stats_subfieldsHc_TLEips_Control_space-T1w_desc-subfields_atlas-bigbrain_dseg_volumetry.png',
#    log: 'logs/hippunfold_subfields_T1w_dseg_volumetry/group_hippunfold_subregions_volumes_tTest_stats.log'
#    group: 'grp_volumes'
#    conda: '../envs/sklearn.yml'
#    script: '../scripts/group_hippunfold_subregions_volumes_tTest_stats.py'


