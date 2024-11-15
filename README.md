# meal_cgm

Preprocessing of CGM data, run the following script one by one (`./cgm/preprocessing`). 

```
preprocess.py
preprocess.r
qc_plot.m
```
Feature extraction, run the `./cgm/feature_extraction/manual_feature.py` and plot by `/cgm/visualization/figures_cgm_features.py`

Analysis, run `run_analysis.r`

## Dependencies
Please refer to the [supplementary files]() for the R and Python environment. The scripts have been tested on Mac OS.

## Installation
Please clone the code and run.

## Demo
We provide a demo on CGM processing and feature extraction. The CGM data can be found [here](). The user can run `cgm/preprocessing/preprocess.r`, `cgm/feature_extraction/manual_feature.py`, and `cgm/visualization/figures_cgm_features.py` in order. The final figures can be compared with Figure 1 and S3. Required time for downloading and running should be brief.
