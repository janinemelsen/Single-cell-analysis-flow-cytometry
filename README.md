# Single-cell-analysis-flow-cytometry
To reveal the cellular heterogeneity within flow cytometric data, clustering analysis and pseudotime analysis can be performed. Visualization can performed by dimensionality reduction. Here we demonstrate the different steps required to perform these analyses.

### 1. Data preprocessing (compensation, export population, transformation and normalization)
[R script](scripts/CSV_to_transformed_normalized_FCS_git.R)


### 2. Clustering, dimensionality reduction, pseudotime
[R script](scripts/clustering_dimensionalityreduction_pseudotime_git.R)


As an example we applied our workflow on fcs files downloaded from the flowRepository (dataset FR-FCM_ZYQ9). Data was compensated in Kaluza, and the live single CD3+ T cells were exported as csv, and further prepocessed by use of the [script](scripts/CSV_to_transformed_normalized_FCS_git.R). We uploaded [the transformed and normalized fcs files](transformed_normalized_CD3/).
Those fcs files were used as input for HSNE-based Gaussian Mean Shift Clustering in Cytosplore. We clustered the CD4 T cells, and exported the clusters as [fcs files](HSNE_clusters_CD4/)

By use of [this script](scripts/clustering_dimensionalityreduction_pseudotime_git.R) we performed alternative dimensionality reduction methods (diffusion map, umap), alternative clustering methods (flowsom, phenograph) and inferred a cellular trajectory (slingshot). The results of this script are saved in [this data frame](df.csv)

The results can be reproduced by loading the [clustered fcs files](HSNE_clusters_CD4/) and use of this [R markdown file](markdown_files/markdown_clustering_dimensionalityreduction_pseudotime.Rmd)

The figures can be reproduced by loading the [data frame](df.csv) and use of this [R markdown file](markdown_files/markdown_visualization.Rmd) 

To view the knitted markdown files click [here](markdown/)

### Downloading the files
To  download all the files go the green button at the right top of the page, and click Download ZIP. If you would like to use our (clustered) fcs files, please note the location where you store the files. To run the scripts, you need to provide the directory where the files are stored.







