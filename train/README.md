# MinNet Training Process
Here is a demonstration on training a MinNet using the human bone marrow mononuclear cells (BMMC) data. Trained model is able to be applied to any BMMC or PBMC data integration.

## Pipeline

### Preparing training dataset
* The training dataset we used to train model is from [GSE194122](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194122) processed H5AD file.
For results in manuscript figures, we kept only batches **s1d1, s1s3, s2d1, s2d4, s2d5, s3d3, s3d6, and s3d10**.
* Save the filtered H5AD file as either `10X/data/10X_Multiome_GEX|ATAC.filt.h5ad` or `cite/data/Cite_GEX|ADT.filt.h5ad`

### Training data processing
* Run the `1_training_data_processing.ipynb` in Jupyter Notebook.
* It will filtering features and cells needed in our previous training.
* It will generate the KNN graph for contrastive loss.
* Save the processed H5AD as either `10X/data/10X_Multiome_GEX|ATAC.training.h5ad` or `cite/data/Cite_GEX|ADT.training.h5ad`

### Training
* Run the `2_training_main.ipynb` in Jupyter Notebook.
* It will perform the training process with default hyper-parameter (the parameter we used).
* Trained weights are saved in `ckpts/`
* During training, histogram of distances between contrastive pairs are showed every 5 epochs.
* Details are written in the jupyter notebook.



## Hyper-paramters
The training process is quite standardized for BMMC situations. But in case users want to train their own model on a specific target tissue, hyper-parameters to be tuned are listed below including **Margin, weights of contrastive loss, and learning rate**.

### Margin


