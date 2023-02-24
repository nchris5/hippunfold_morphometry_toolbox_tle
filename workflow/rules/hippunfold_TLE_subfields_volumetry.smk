#Volumetry
rule volumetry_T1w_subfields:
    input:
        dseg_T1w_subfields = config['dseg_subfields_bigbrain_spaceT1w'],
    params:
        voxel_shape = config['voxel_shape'],
    output:
        volumetry_T1w_subfields = 'work/hippunfold_volumetry/sub-{subject}/sub-{subject}_hemi-{hemi}_space-T1w_desc-subfields_atlas-bigbrain_volumes.tsv',
    log: 'logs/hippunfold_volumetry/sub-{subject}_hemi-{hemi}_volumetry_T1w_subfields_dseg.log'
    group: 'VolmtyLR'
    script: '../scripts/volumetry_subfields.py'

#Group individual subjects hippocampal subfields into a single tsv file
rule group_volumetry_T1w_subfields:
    input:
        L_volumetry_T1w_subfields_tsv = expand(rules.volumetry_T1w_subfields.output, subject=subjects, hemi='L'),
        R_volumetry_T1w_subfields_tsv = expand(rules.volumetry_T1w_subfields.output, subject=subjects, hemi='R'),
    params:
        demographics_tsv = config['demographics_tsv'],
    output:
        group_volumetry_T1w_subfields_tsv_dir = directory('results/hippunfold_volumetry/group_volumetry/'),
    log: 'logs/hippunfold_volumetry/group_volumetry_T1w_subfields_dseg.log'
    group: 'grpVolmty'
    script: '../scripts/group_volumetry_subfields.py'

rule tTest_group_volumetry_T1w_subfields:
    input:
        group_volumetry_T1w_subfields_tsv_dir = rules.group_volumetry_T1w_subfields.output,
    params: #Used to get data for each subgroup and subfield within the py script
        group_names = subgroups,
        group_sf_volumetry_names = ['Subiculum_volume', 'CA1_volume', 'CA2_volume', 'CA3_volume', 'CA4_volume', 'DG_volume'],
    output:
        stats_group_volumetry_T1w_subfields_tsv_dir = directory('results/hippunfold_volumetry/group_volumetry_stats/'),
    log: 'logs/hippunfold_volumetry/stats_group_volumetry_T1w_subfields_dseg.log'
    group: 'grpVolmty'
    script: '../scripts/stats_group_volumetry_subfields.py'


