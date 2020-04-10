# Single-cell-analysis-flow-cytometry
To reveal the cellular heterogeneity within flow cytometric data, clustering analysis and pseudotime analysis can be performed. Visualization can performed by dimensionality reduction. Here we demonstrate the different steps required to perform these analyses.

1. Data preprocessing (compensations, export population, transformation and normalization.

--> CSV_to_transformed_normalized_fcs R script

2. Clustering, dimensionality reduction, pseudotime

--> clustering_dimensionalityreduction_pseudotime R script
--> clustering_dimensionalityreduction_pseudotime R markdown

As an example we uploaded the transformed and normalized fcs files from dataset FR-FCM-ZYQ9, as published at the flowRepository. 

Data was clustered by HSNE-based Gaussean Mean shift (GMS) clustering in Cytosplore. Each cluster was exported as fcs file and uploaded.
The data frame df, includes all tranformed and normalized expression values, clustering assignments by flowsom, phenograph and HSNE-based GMS, reduced dimensions (diffusion map, umap), and pseudotime values.

By applying the R markdown document 'visualization' figures can be reproduced.





