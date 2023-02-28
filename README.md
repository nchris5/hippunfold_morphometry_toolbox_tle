# hippunfold_morphometry_toolbox_tle
Morphometic and Volumetric Analyses using hippunfold in patients with temporal lobe epilepsy; Primarily investigating CA1 inferior hippocampal digitations


1) Calculation of displacement metric from native hippocampal surface to smoothed native surface; displacment metric is then propagated to unfolded volumetric space for subsequent analyses.

2) Curve generation of anterior-posterior inferior CA1 hippocampal digitations in unfolded volumetric space by averaging CA1 displacement at each level of the unfolded anterior-posterior dimensiom in the hippocampal body's outer surface.
  -  This provides a single-dimensional (AP) signal permitting quantification of the extent of hippocampal digitation; fourier analyses to quantify the frequency and amplitude of digitations are subsequently performed.
  -  Similar analyses conducted for any/all of the following additional hippunfold metrics provided (Curvature, Gyrification, Thickness, Surfarea).

3) Volumetry of the whole hippocampus and individual subfields (Subiculum, CA1, CA2, CA3, CA4, DG).

4) All previous metrics are grouped into the following possible subgroups for subsequent statistical analysis.
  -  ipsilateral_TLE,  contralateral_TLE,  bilateral_TLE, allAffected_TLE(ipsilateral+bilateral),  all_TLE(ipsilateral+contralateral+bilateral), controls

---INDEV---

5) Quantification and group statistical analyses of CA1 inferior hippocampal digitations through fourier analyses
-  A) Frequency (number of digitations)
-  B) Amplitude (size of digitations)

---INDEV---


Installation Instructions:
1) Clone rep with
- git clone https://github.com/nchris5/hippunfold_morphometry_toolbox_tle
2) Poetry install and activate venv
- cd hippunfold_morphometry_toolbox_tle
- poetry install
- poetry shell
