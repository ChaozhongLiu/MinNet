# Model-based smoothing and cis-regulatory element inferring
Here is the codes for reproducing manuscript Figure 4 and the supplementary figures related to it.


## Data availability
The single-cell data used in figure 4 results is the 10X Multiome PBMC dataset from 10X Genomics (same in figure 2).  
The pcHi-C evidence is download from https://ars.els-cdn.com/content/image/1-s2.0-S0092867416313228-mmc4.zip and https://osf.io/e594p/. After downloading, save the data in a folder called `pcHiC` here.

## Procedure
* Run `data/Pre-processing_pcHiC.ipynb` to prepare everything needed for pcHiC evidence.
* Run `SiaNN_PBMC.ipynb` to generate embedding space needed below.
* Run `data/Pre-processing_smoothing.ipynb` to prepare PBMC data for smoothing.
* Run `data/PBMC_smooth.ipynb` to generate smoothed data with user-provided number of neighbors.

## Evaluation
* Run `PBMC_smooth_summary.ipynb` to make sparsity summary, correlation between genes and 2kb nearby peaks summary.
* Run `PBMC_smooth_pcHi-C.ipynb` to validate correlation between genes and 150kb nearby peaks with pcHi-C evidence.
* Run `Plots.R` for summary plots, `genome_track.R` for example genome track visualization.


