# Model studay case: COVID-19 singe-cell dataset
Here is the codes for reproducing manuscript Figure 5 and the supplementary figures related to it.


## Data availability
The COVID-19 dataset is from the paper [Multi-omic profiling reveals widespread dysregulation of innate immunity and hematopoiesis in COVID-19](https://rupress.org/jem/article/218/8/e20210582/212379/Multi-omic-profiling-reveals-widespread). 
The pre-processed scRNA-seq data can be downloaded at their [GitHub](https://hosted-matrices-prod.s3-us-west-2.amazonaws.com/Single_cell_atlas_of_peripheral_immune_response_to_SARS_CoV_2_infection-25/Single_cell_atlas_of_peripheral_immune_response_to_SARS_CoV_2_infection.h5ad) in H5AD format





## Procedure
* Run `data/Pre-processing_pcHiC.ipynb` to prepare everything needed for pcHiC evidence.
* Run `SiaNN_PBMC.ipynb` to generate embedding space needed below.
* Run `data/Pre-processing_smoothing.ipynb` to prepare PBMC data for smoothing.
* Run `data/PBMC_smooth.ipynb` to generate smoothed data with user-provided number of neighbors.

## Evaluation
* Run `PBMC_smooth_summary.ipynb` to make sparsity summary, correlation between genes and 2kb nearby peaks summary.
* Run `PBMC_smooth_pcHi-C.ipynb` to validate correlation between genes and 150kb nearby peaks with pcHi-C evidence.
* Run `Plots.R` for summary plots, `genome_track.R` for example genome track visualization.


