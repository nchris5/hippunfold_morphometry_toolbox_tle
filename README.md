# hippunfold_morphometry_toolbox_tle
Morphometic and Volumetric Analyses using hippunfold in patients with temporal lobe epilepsy; Primarily investigating CA1 digitations


1) Calculation of native unfolded hippocampal outer surface displacement relative to smoothed surface

2) Curve generation of anterior-posterior inferior CA1 hippocampal digitations in unfolded volumetric space by averaging CA1 displacement at each level of the anterior-posterior dimensiom in the hippocampal body's outer surface.
  -  This provides a single-dimensional (AP) signal permitting quantification of the extent of hippocampal digitation; fourier analyses are subsequently performed

3) Similar analyses are performed for hippocampal curavture and gyrification, in addition to subfield volumetry

4) All previous metrics are grouped into the following subgroups of interest for subsequent statistical analysis
  -  ipsilateral_TLE,  contralateral_TLE,  bilateral_TLE,  allAffected_TLE(ipsilateral+bilateral),  all_TLE(ipsilateral+contralateral+bilateral),  controls,
  -  These subgrouped metrics are provided in the results directory as tsv files with each row representing a participant and their respective metric values



Dependencies: (see https://hippunfold.readthedocs.io/en/latest/getting_started/installation.html)
  On the cbs_server, best solution is to install and activate the hippunfold virtual environment as follows:
1) Clone hippunfold with
git clone https://github.com/khanlab/hippunfold
2) Install and activate hippunfold venv
cd hippunfold
poetry install
poetry shell
