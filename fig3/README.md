# MinNet Batch effect removal
Here is the codes for reproducing manuscript Figure 3 and the supplementary figures related to it.

## Directory Structure
This directory is organized by algorithms.
* `data/data_preprocessing.ipynb` is for testing data generation.
* Algorithms' results are all saved in `dataset/results/raw/`
* After running all algorithms, summary evaluation was done using `evaluation.ipynb`. Results were saved in `dataset/results/`
* To save space, all datasets and raw results were deleted.
* Below are specific explanations on how testing cases were setted and evaluated.

**FYI**: during development, our algorithm used the temperory name 'SiaNN', that's why in all codes SiaNN is seen.

## Reproduction Pipeline
All the testing evaluations follows the same pipeline:
* Data pre-processing
* Run all algorithms in the order of MinNet - Seurat - others
* Run `dataset/evaluation.ipynb` to get the summary statistics in `dataset/results/`
* Finally, plots were drawn by running batch_removal_plots.R using data in `dataset/results/`.

## Testing scenario generation
### Situation to integrate scRNA-seq and scATAC-seq with the same batches
This is the testing case 1 in fig2 dataset, so no extra implementation is needed. We just took the results during evaluation and plotting.

### Situation to integrate scRNA-seq and scATAC-seq from different batches
The 2 cases simulate the real problem when researchers get their single-cell from different samples, or even different resources. The ability of an integration algotrithm to distinguish between batch effect and actual biological variance is essential here.

Case 1: scRNA-seq data from NeurIPS `s3d7` and scATAC-seq from `s4d1`
Case 2: scRNA-seq data from NeurIPS `s4d1` and scATAC-seq from `s3d7`


## Algorithm implementation
1. MinNet: Run the SiaNN.ipynb in `dataset/SiaNN/`
2. Seurat v3: Run h5ad2rds.R to transfer from python to R. Then run seurat.R followed by after_seurat.ipynb for evaluation metrics.
3. bindSC: Run bindSC.R, followed by after_bindSC.ipynb.
4. GLUE: prepare H5AD and graph for GLUE by runing glue_preprocessing.ipynb. Then glue.py 8 times using glue.sh due to the randomness in algorithm.
5. All versions of Liger: run liger.R and after_liger.py 8 times by the shell scripts in each folder.




