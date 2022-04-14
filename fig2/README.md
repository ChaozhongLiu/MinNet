# MinNet Benchmarking
Here is the codes for reproducing manuscript Figure 2 and the supplementary figures related to it.

## Directory Structure
This directory is organized by datasets.
* Within each dataset's folder, codes are organized by algorithms, e.g., `dataset/seurat/`
* Algorithms' results are all saved in `dataset/seurat/results/raw/`
* After running all algorithms, summary evaluation was done using `evaluation.ipynb`. Results were saved in `dataset/seurat/results/`
* To save space, all datasets and raw results were deleted.
* Below are specific instructions on how to get the data and get the raw results.

FYI: during development, our algorithm used the temperory name 'SiaNN', that's why in all codes SiaNN is seen.

## Reproduction Pipeline
All the testing evaluations follows the same pipeline:
* Data downloading and pre-processing
* Run all algorithms in the order of MinNet - Seurat - others
* Run `dataset/evaluation.ipynb` to get the summary statistics in `dataset/seurat/results/`
* Finally, plots were drawn by running benchmark_plots.R using data in `dataset/seurat/results/`.

## Datasets Downloading
### 10X Multiome datasets and Cite-seq datasets of BMMC
Datasets names as BMMC_data, BMMC_test, Cite_data, and Cite_data_2 are from the [NeurIPS 2021 Competition](https://openproblems.bio/neurips_2021/).
The pre-processed H5AD file can be downloaded at [GSE194122](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194122).
The batch spliting is:
* BMMC_data: `s1d2 and s3d7`
* BMMC_test: `s4d1, s4d8, s4d9`
* Cite_data: `s1d2 and s3d7`
* Cite_data_2: `s4d1, s4d8, s4d9`


### 10X Multiome PBMC dataset
The dataset can be downloaded at 10X Genomics website [here](https://support.10xgenomics.com/single-cell-multiome-atac-gex/datasets/1.0.0/pbmc_granulocyte_sorted_10k)

### 10X Multiome Human brain dataset
The dataset can be downloaded at 10X Genomics website [here](https://www.10xgenomics.com/resources/datasets/frozen-human-healthy-brain-tissue-3-k-1-standard-2-0-0)


## Algorithm implementation
1. MinNet: Run the SiaNN.ipynb in `dataset/SiaNN/`
2. Seurat v3: Run h5ad2rds.R to transfer from python to R. Then run seurat.R followed by after_seurat.ipynb for evaluation metrics.
3. bindSC: Run bindSC.R, followed by after_bindSC.ipynb.
4. GLUE: prepare H5AD and graph for GLUE by runing glue_preprocessing.ipynb. Then glue.py 8 times using glue.sh due to the randomness in algorithm.
5. All versions of Liger: run liger.R and after_liger.py 8 times by the shell scripts in each folder.




