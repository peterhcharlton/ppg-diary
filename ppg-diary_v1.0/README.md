# PPG Diary Scripts: v.1

This version of the toolbox contains the scripts used in the following publication:

Charlton P.H. *et al.* **Acquiring wearable photoplethysmography data in daily life: the PPG Diary Pilot Study**, in *Engineering Proceedings*, [Under Review]

The scripts are provided in Matlab &reg; format.

## Summary of Publication

This pilot study investigated the feasibility of acquiring photoplethysmography (PPG) data in daily living. Several lessons were learnt about how to acquire PPG data in future studies.

## Replicating this Publication

The work presented in this case study can be replicated as follows:

*   Download the *ppg_diary_pilot_conv_data.mat* data file from [here](https://doi.org/10.5281/zenodo.3268500).
*   Download [Version 1](https://github.com/peterhcharlton/ppg-diary/tree/master/ppg-diary_v1.0) of the scripts.
*   Specify the path of the *ppg_diary_pilot_conv_data.mat* data file in line 111 of the *ppg_diary_pilot1_analysis.m* script.
*   Analyse the data using the *[ppg_diary_pilot1_analysis.m](https://raw.githubusercontent.com/peterhcharlton/ppg-diary/master/ppg-diary_v1.0/ppg_diary_pilot1_analysis.m)* script. This script calls the *[PulseAnalyse.m](https://raw.githubusercontent.com/peterhcharlton/ppg-diary/master/ppg-diary_v1.0/PulseAnalyse.m)* script, so you will need this script too (v.1.2 beta).

