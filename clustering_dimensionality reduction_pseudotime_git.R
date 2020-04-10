###############################
## INSTALL REQUIRED PACKAGES ##
###############################
install.packages("reshape2")
install.packages("ggplot2")
install.packages("uwot")
install.packages("ggrepel")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("scales")
install.packages("reshape2")
install.packages("RColorBrewer")
install.packages("devtools")
install.packages("BiocManager")

BiocManager::install("FlowSOM")
BiocManager::install("slingshot")
BiocManager::install("flowCore")
BiocManager::install("SingleCellExperiment")

library(devtools)
devtools::install_github("JinmiaoChenLab/cytofkit2")
devtools::install_github('flying-sheep/knn.covertree')
devtools::install_github('theislab/destiny')

###############################
### LOAD REQUIRED PACKAGES ####
###############################

library(flowCore)
library(FlowSOM)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(scales)
library(reshape2)
library(RColorBrewer)
library(destiny)
library(uwot)
library(slingshot)
library(cytofkit2)
library(ggrepel)

############################
######## LOAD DATA #########
############################

### Load the (transformed, normalized, unclustered) FCS files from the 'CSV_to_transformed_normalized_FCS' script
### Or load fcs files which where clustered in Cytosplore (in this case each fcs file is 1 cluster)

## Provide the directory of the fcs files
dirFCS = "path_to_your_folder_with_fcs_files" 

## Optional: when loading clustered fcs files from cytosplore, provide the directory of the text file 'CSPLR_ST.txt'. Cytosplore exports this file upon running the HSNE. This file contains the decoding of the sample numbers.
pathST <- "I:/Laboratoria/immlab/Janine/flowpaper/JI/rebuttal/newfigures/CD4_level2/clusters/CSPLR_ST.txt"

## Defining a function to read multiple fcs files from a directory 'dir' into a single data.frame:
# NB: The column in the output named 'fileName' tracks the original file where each cell came from.
# Optionally perform remapping of column 'CSPLR_ST' holding cytosplore sample numbers to actual names:
read.flowdat <- function(dir,path_CSPLR_ST=""){
  # Read:
  filepaths <- list.files(path=dir,pattern = ".fcs", full.names=TRUE)
  flowset <- read.flowSet(files=filepaths, transformation=FALSE, truncate_max_range = FALSE)
  # Transform to data frame:
  x <- as.data.frame(exprs(as(flowset,'flowFrame')),stringsAsFactors=FALSE)
  # Map column 'Original' to filename (in this case holding clusters of HSNE):
  filenames <- gsub("[.fcs]","",list.files(path=dir,pattern = ".fcs", full.names=FALSE))
  names(filenames) <- sort(unique(x$Original))
  x$fileName <- filenames[as.character(x$Original)]
  # Remove column 'Original':
  x <- x[,-which(colnames(x)=="Original")]
  # Optionally remap Cytosplore sample tags to original filename:
  if(file.exists(path_CSPLR_ST)){
    # Read:
    sampID <- gsub(".fcs","",basename(sapply(strsplit(readLines(path_CSPLR_ST),": "),function(x) x[1])))
    names(sampID) <- sapply(strsplit(readLines(path_CSPLR_ST),": "),function(x) x[2])
    x$sampleID <- sampID[as.character(x$CSPLR_ST)]
  }
  return(x)
}

## Read fcs files
# In our example we will read the data which were clustered in Cytosplore (each fcs file is 1 cluster)
df <- read.flowdat(dir=dirFCS,path_CSPLR_ST = pathST)
# Optional: Set columname 'fileName' to clusters_HSNE:
colnames(df)[which(colnames(df)=="fileName")] <- "clusters_HSNE"

## In our example we will start with the 275856 CD4 T cells

############################
######## CLUSTERING ########
############################

### We will discuss 3 clustering methods:
### A) HSNE-based Gaussian Mean Shift clustering 
### B) flowSOM
### C) Phenograph

## Option A: HSNE-based Gaussian Mean Shift clustering
## With the software cytosplore HSNE-based GMS clustering can be performed. The clustering results with the expression values can be exported as fcs and loaded as described above

## Option B: generate clusters by FlowSOM
#check colnames, to determine which columns you need for the cluster calculation
colnames(df)
#run flowsom
flowsom <- FlowSOM(input = dirFCS, 
                transform = FALSE,
                scale = FALSE,
                colsToUse = c(7:9, 11, 13:16,18,19), #provide the columns for the clustering
                nClus = 14, #we choose 14, since we also generated 14 clusters by HSNE
                seed = 100)

# Get metaclustering per cell
clusters_flowsom <- as.factor(flowsom$FlowSOM$map$mapping[,1])
levels(clusters_flowsom) <- flowsom$metaclustering

#add flowsom clusters to dataframe
df <- cbind(df, clusters_flowsom)

## Option C: generate clusters by Phenograph (based on Louvain clustering)
# select the columns for the clustering calculation (usually the numbers are the same as used for the flowsom calculation)
#the higher the K nearest neighbours, the lower the number of clusters
phenograph <- Rphenograph(df[,c(7:9,11, 13:16,18,19)], k=50)
clusters_phenograph <- as.factor(phenograph$membership)

#add phenograph clusters to expression data frame
df <- cbind(df, clusters_phenograph)

###############################
######## VISUALIZATION ########
###############################
### We will discuss 3 visualization methods:
### A) HSNE
### B) diffusion map
### C) UMAP

## option A: HSNE
# HSNE can only be visualized in Cytosplore. Import the FCS files from the 'CSV_to_transformed_normalized_FCS' script to cytosplore to perform HSNE analysis

## Option B: diffusion map (for our example CD4 T cell dataset this will take approximately 2hours)
# reduce the K, if computational load is too high [it takes approximately 2 hours for the example dataset of 275856 cells]
dm <- DiffusionMap(df, vars = c("CD95", "CD8", "CD27", "CCR7", "CD45RA", "CD49B", "CD69", "CD103", "CD3", "CD4"), k=1000, suppress_dpt = TRUE, verbose=TRUE)

# add the diffusion components to the expression data frame (either all by dm@eigenvectors, or a selection by dm$DC1, dm$DC2, etc.)
df <- cbind(df, DC1=dm$DC1, DC2=dm$DC2, DC3=dm$DC3)

# Vizualize results:
# With the viz.dm function the two components of the diffusion map can be visualized with the expression levels of a particular parameter applied as color scale
# dat: the data frame
# param.name: the name of the parameter which you would like to plot as color scale
# if no limits are provided: for better visualization, the most extreme low and high percentile of the color values are normalized.
# if limits are provided: the color scale will be adjusted according to the provided limits (recommended for parameters which are expressed on all cells). To apply the same colorscale as in Cytosplore, provide the limits as set in Cytosplore.

viz.dm <- function(dat,dr, param.name,limits=NULL){
  ColVal <- dat[,param.name]
  if(is.null(limits)){
    Lim <- quantile(ColVal,probs=seq(0,1,0.01))[c(2,100)]
    p <- ggplot(dat, aes(x = DC1, y =DC2)) +geom_point(aes(color = ColVal), size=0.1)+theme_classic()+scale_color_distiller(name=param.name, palette = "RdYlBu", limits=Lim, oob=squish)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+ggtitle(param.name)
  } else {
    p <- ggplot(dat, aes(x = DC1, y = DC2)) +geom_point(aes(color = ColVal), size=0.1)+theme_classic()+scale_color_distiller(name=param.name, palette = "RdYlBu", limits=c(limits[1],limits[2]), oob=squish)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+ggtitle(param.name)
  }
  p
}

viz.dm(dat=df,param.name='CD27', limits=c(-0.02,3.17))
viz.dm(dat=df, param.name="CD45RA")

#visualize and label clusters on diffusion map
#first we need to determine the positions of the cluster labels, based on the DC coordinates
label_HSNE_dm <- df%>%group_by(clusters_HSNE)%>%select(DC1, DC2)%>%summarize_all(mean)
label_flowsom_dm <- df%>%group_by(clusters_flowsom)%>%select(DC1, DC2)%>%summarize_all(mean)
label_pheno_dm <- df%>%group_by(clusters_phenograph)%>%select(DC1, DC2)%>%summarize_all(mean)

ggplot(df, aes(x=DC1, y=DC2, color=clusters_HSNE))+geom_point(size=0.1)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+geom_label_repel(aes(label=clusters_HSNE), data=label_HSNE_dm)+guides(colour=FALSE)
ggplot(df, aes(x=DC1, y=DC2, color=as.factor(clusters_flowsom)))+geom_point(size=0.1)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+geom_label_repel(aes(label=clusters_flowsom), data=label_flowsom_dm)+guides(colour=FALSE)
ggplot(df, aes(x=DC1, y=DC2, color=as.factor(clusters_phenograph)))+geom_point(size=0.1)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+geom_label_repel(aes(label=clusters_phenograph), data=label_pheno_dm)+guides(colour=FALSE)


## Option C: UMAP
# select the columns for the UMAP calculation
# check different n_neighbours (controls how UMAP balances local versus global structure in the data) for your UMAP plot
# check min_dist (controls how tightly UMAP is allowed to pack points together, low values=clumpier embeddings) for your UMAP plot
umap <- umap(df[,c(7:9, 11, 13:16, 18, 19)], n_neighbors = 30, min_dist=0.001, verbose=TRUE)
umap<- as.data.frame(umap)
colnames(umap) <- c('umap_1', 'umap_2')
df <- cbind(df,umap)


# Vizualize results:
# With the viz.umap function the two components of the umap can be visualized with the expression levels of a particular parameter applied as color scale
# dat: the data frame
# param.name: the name of the parameter which you would like to plot as color scale
# if no limits are provided: for better visualization, the most extreme low and high percentile of the color values are normalized.
# if limits are provided: the color scale will be adjusted according to the provided limits (recommended for parameters which are expressed on all cells). To apply the same colorscale as in Cytosplore, provide the limits as set in Cytosplore.
viz.umap <- function(dat,param.name,limits=NULL){
  ColVal <- dat[,param.name]
  if(is.null(limits)){
    Lim <- quantile(ColVal,probs=seq(0,1,0.01))[c(2,100)]
    p <- ggplot(dat, aes(x = umap_1, y =umap_2)) +geom_point(aes(color = ColVal), size=0.1)+theme_classic()+scale_color_distiller(name=param.name, palette = "RdYlBu", limits=Lim, oob=squish)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+ggtitle(param.name)
  } else {
    p <- ggplot(dat, aes(x = umap_1, y = umap_2)) +geom_point(aes(color = ColVal), size=0.1)+theme_classic()+scale_color_distiller(name=param.name, palette = "RdYlBu", limits=c(limits[1],limits[2]), oob=squish)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+ggtitle(param.name)
  }
  p
}

viz.umap(dat=df,param.name='CD27', limits=c(-0.02,3.17))
viz.umap(dat=df, param.name="CD45RA")

#visualize and label clusters on umap
label_HSNE_umap <- df%>%group_by(clusters_HSNE)%>%select(umap_1, umap_2)%>%summarize_all(mean)
label_flowsom_umap <- df%>%group_by(clusters_flowsom)%>%select(umap_1, umap_2)%>%summarize_all(mean)
label_pheno_umap <- df%>%group_by(clusters_phenograph)%>%select(umap_1, umap_2)%>%summarize_all(mean)

ggplot(df, aes(x=umap_1, y=umap_2, color=as.factor(clusters_HSNE)))+geom_point(size=0.1)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+geom_label_repel(aes(label=clusters_HSNE), data=label_HSNE_umap)+guides(colour=FALSE)
ggplot(df, aes(x=umap_1, y=umap_2, color=as.factor(clusters_flowsom)))+geom_point(size=0.1)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+geom_label_repel(aes(label=clusters_flowsom), data=label_flowsom_umap)+guides(colour=FALSE)
ggplot(df, aes(x=umap_1, y=umap_2, color=as.factor(clusters_phenograph)))+geom_point(size=0.1)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+geom_label_repel(aes(label=clusters_phenograph), data=label_pheno_umap)+guides(colour=FALSE)

#############################
######## FREQUENCIES ########
#############################
# As an example, we show the cluster frequencies of the HSNE clusters. To check the other cluster frequencies, adjust df$clusters_HSNE
# check frequency of each cluster per sample
counts <- as.data.frame.matrix(table(df$sampleID,df$clusters_HSNE))

#calculate percentages of a certain sample present in cluster
counts_percofsample = counts/rowSums(counts)*100
counts_percofsample$total <- rowSums(counts)
counts_percofsample$sampleID <- row.names(counts_percofsample)
counts_percofsample <- melt(counts_percofsample,  id.vars=c('sampleID', 'total'), variable.name='cluster', value.name='frequency')

ggplot(counts_percofsample, aes(x=reorder(cluster, frequency, FUN=median), y=frequency))+ geom_boxplot(position='dodge2', outlier.shape=NA)+geom_jitter(aes(colour=sampleID), size=2)+xlab('cluster')
ggplot(counts_percofsample, aes(x=cluster, y=frequency, fill=cluster))+geom_bar(stat='identity')+facet_wrap(~sampleID)+theme(axis.text.x=element_text(angle=30, vjust=.9, hjust=1))


############################
#### COMPUTE SLINGSHOT #####
############################

## To calculate pseudotime we will use Slingshot
## Slingshot requires clusters as input for the lineage identification. In this example, we use the HSNE based clusters.
## Too many clusters will lead to artificial lineages, therefore we merge the CD49b- and CD49b+ subset of each naive and memory population, and relabel them
df$merged_HSNE <-gsub("CD4-1$", "Naive", df$clusters_HSNE)
df$merged_HSNE <-gsub("CD4-2$", "Naive", df$merged_HSNE)
df$merged_HSNE <-gsub("CD4-3$", "CD27+ CM", df$merged_HSNE)
df$merged_HSNE <-gsub("CD4-4$", "CD27+ CM", df$merged_HSNE)
df$merged_HSNE <-gsub("CD4-5$", "CD27+ EM", df$merged_HSNE)
df$merged_HSNE <-gsub("CD4-6$", "CD27+ EM", df$merged_HSNE)
df$merged_HSNE <-gsub("CD4-7$", "CD27- EM", df$merged_HSNE)
df$merged_HSNE <-gsub("CD4-8$", "CD27- EM", df$merged_HSNE)
df$merged_HSNE <-gsub("CD4-9$", "EMRA", df$merged_HSNE)
df$merged_HSNE <-gsub("CD4-10$", "EMRA", df$merged_HSNE)
df$merged_HSNE <-gsub("CD4-11$", "CD8dim EMRA", df$merged_HSNE)
df$merged_HSNE <-gsub("CD4-12$", "CD8dim EM", df$merged_HSNE)
df$merged_HSNE <-gsub("CD4-13$", "CD69+CD103- EM", df$merged_HSNE)
df$merged_HSNE <-gsub("CD4-14$", "CD103+ EM", df$merged_HSNE)


#to visualize diffusionmap with merged clusters
label_merged_HSNE_dm <- df%>%group_by(merged_HSNE)%>%select(DC1, DC2)%>%summarize_all(mean)
ggplot(df, aes(x=DC1, y=DC2, color=merged_HSNE))+geom_point(size=0.1)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+geom_label_repel(aes(label=merged_HSNE), data=label_merged_HSNE_dm)+guides(colour=FALSE)

##Create a Slingshot object and add cell expression and cluster information to this object
#create slingshot object
slingshot_object <- SingleCellExperiment(assays = List(norm = as.matrix(t(df))))

#add expression data to slingshot and select markers to use for the pseudotime calculation
reducedDims(slingshot_object) <- SimpleList(expressiondata = as.matrix(df%>%select(c("CD95", "CD8", "CD27", "CCR7", "CD45RA", "CD49b", "CD69", "CD103", "CD3", "CD4"))))

#Add clusters to slingshot object. Either the HSNE-based clusters, the Phenograph clusters, or flowsomclusters.
#In this example we demonstrate the use of the merged HSNE-based clusters.
colData(slingshot_object)$clusters <- df$merged_HSNE

#calculate lineages (optional: appoint starting cluster)
#In this example we appoint the naive T cells, as starting cluster.
lin <- getLineages(reducedDims(slingshot_object)$expressiondata,colData(slingshot_object)$clusters, start.clus='Naive')

#construct smooth curves and calculate pseudotime [This will take approximately 72 hours for this example dataset]
#to reduce computational time, curves can be approximated by a fixed number of points, for instance 100.
curve <- getCurves(lin, approx_points=FALSE)

#add pseudotime values to the data frame
df <- cbind(df, as.data.frame(slingPseudotime(curve), row.names=FALSE))

#generate table with diffusion map coordinates, pseudotimevalues and cluster assignment (either the HSNE-based clusters or other clusters)
#In this example we use the HSNE-based clusters
#If you would like to use other clusters, adjust 'merged_HSNE' to clusters_flowsom', for instance
pseudotimevalues <- df%>%select(c(DC1, DC2,clusterID=merged_HSNE, curve1, curve2, curve3))

#reshape table 
pseudotimevalues <- melt(pseudotimevalues, id.vars=c('clusterID', 'DC1', 'DC2'), variable.name='lineage', value.name = 'pseudotime')

#rename curve to lineage
pseudotimevalues$lineage <- gsub('curve','lineage', pseudotimevalues$lineage)

#exclude cells with NA pseudotimevalues (those cells are not present in all lineages and have NA values for the lineages in which they are absent)
pseudotimevaluesexclNA <- pseudotimevalues%>%filter(pseudotime!='NA')

#generate colorpalette
colors <- colorRampPalette(rev(brewer.pal(11, 'Spectral'))[-6])(100)

#ggplot of each lineage colored by either pseudotime
ggplot(pseudotimevaluesexclNA%>%arrange(pseudotime), aes(x=DC1, y=DC2))+geom_point(aes(color=pseudotime),size=0.1, alpha=0.3) +facet_wrap(~lineage)+scale_color_gradientn(colours=colors)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())

#ggplot of each lineage colored by clusters (in our example, clusterID)
ggplot(pseudotimevaluesexclNA%>%arrange(pseudotime), aes(x=DC1, y=DC2))+geom_point(aes(color=clusterID),size =0.1)+facet_wrap(~lineage)+theme_bw()

#jitterplot of pseudotimevalues per lineage
ggplot(pseudotimevaluesexclNA%>%arrange(pseudotime), aes(x=pseudotime, y=lineage)) + geom_jitter(aes(color=clusterID), size=0.1)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())

### save final df file  ###
write.csv(df, 'df.csv')
