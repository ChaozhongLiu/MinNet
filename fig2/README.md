# MinNet Benchmarking
Here is all the codes used for reproducing manuscript Figure 2 and the supplementary figures related to it.

## Directory Structure
This directory is organized by datasets.
* Within each dataset's folder, codes are organized by algorithms, e.g., `dataset/seurat/`
* Algorithms' results are all saved in `dataset/seurat/results/raw/`
* After running all algorithms, summary evaluation was done using `evaluation.ipynb`. Results were saved in `dataset/seurat/results/`
* To save space, all datasets and raw results were deleted.
* Specific instructions on how to get the data and get the raw results are below.


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
