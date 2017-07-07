# app-AFQ-seg
Segments White Matter Tracts

## This service currently runs on following computing resources

### IU Carbonate

It requires a conda environment created with following packages installed

```
conda create -n app-afq-seg
source activate app-afq-seg
conda config --add channels conda-forge
conda install boto3
conda install numpy
conda install scipy
```

On all computing resources, following packages must be installed and accessible 

```
/N/u/aryaam/Karst/dipy
/N/u/aryaam/Karst/github_repos/nibabel
/N/u/aryaam/Karst/github_repos/pyAFQ
```
