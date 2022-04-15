# Model studay case: COVID-19 singe-cell dataset
Here is the codes for reproducing manuscript Figure 5 and the supplementary figures related to it.



## Data availability
The COVID-19 dataset is from the paper [Multi-omic profiling reveals widespread dysregulation of innate immunity and hematopoiesis in COVID-19](https://rupress.org/jem/article/218/8/e20210582/212379/Multi-omic-profiling-reveals-widespread).  
The pre-processed scRNA-seq data can be downloaded [here](https://hosted-matrices-prod.s3-us-west-2.amazonaws.com/Single_cell_atlas_of_peripheral_immune_response_to_SARS_CoV_2_infection-25/Single_cell_atlas_of_peripheral_immune_response_to_SARS_CoV_2_infection.h5ad) in H5AD format and save in `data/scRNA/`.  
The scATAC-seq data is downloaded from [GSE174072](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174072) and should be saved in `data/scATAC/raw/`.  

The samples we used are listed in `data/scATAC/raw/sample_list.txt`.  
The DEGs needed for analysis are provided in `data/`.



## Procedure
### PBMC model training
* Processing the 10X Multiome PBMC data by running `PBMC_model_train/data/data_processing.ipynb`.
* Run `train_main.ipynb` to train the model with alaready setted parameters.
* In `ckpts/`, find the `best_model_margin2.5_lamb1.0.ckpt` and save as `../utils/SiaNN.PBMC.ckpt`.
* The trained weights already saved in `../utils/SiaNN.PBMC.ckpt` if you want to save time.

### COVID-19 data processing
* For scRNA-seq data, run `scRNA/data_processing.ipynb`.
* For scATAC-seq data, run `scATAC/raw/JEM_ATAC_process_1/2.R`, followed by `scATAC/JEM_ATAC_process_3.ipynb`.
* Finally, run `data/consistence_processing.ipynb` to re-annotate cell types in scRNA-seq and scATAC-seq.
* The final processed data will be saved in `data/JEM_GEX.h5ad and JEM_ATAC.h5ad`.

### Applying MinNet
* Run `SiaNN_BMMCmodel.ipynb` to apply BMMC-trained model and calculate embedding space, label transferring, etc.
* Run `SiaNN_PBMCmodel.ipynb` to apply PBMC-trained model and calculate embedding space, label transferring, etc.


### Cis-regulatory element inferring
Inferring was done per cell type. We chose NK and Monocytes as the target cell types.
* Before calculating spearman correlation, run `Smooth_jem_covid.ipynb` for data smoothing.
* Run `cis_regulatory_inferring.ipynb` to calculate the gene-peak correlation. Summary and mini-bulk data are saved in `cis_infer`.

### Severity-stratified correlation analysis
We chose NK and Monocytes as the target cell types.
* Run `cor_by_severity.ipynb` to calculate the gene-peak correlation per severity group.
* Summary and mini-bulk data are saved in `cor_by_severity`.

## Evaluation
* Run `Plots.R` to plot the UMAP visualization and label transfer accuracy heatmap.
* Run `Plots_2.R` to plot correlation related analysis.
* Run `genome_track.R` to plot example case genome track.


