###############################
## INSTALL REQUIRED PACKAGES ##
###############################

install.packages("BiocManager")
install.packages("dplyr")
BiocManager::install("Biobase")
BiocManager::install("flowCore")
BiocManager::install("flowVS")

###############################
### LOAD REQUIRED PACKAGES ####
###############################

library(flowCore)
library(Biobase)
library(dplyr)
library(flowVS)

###############################
##### LOAD PREPARED DATA ######
###############################

## Save the the compensated events in the gate of interest per individual sample as a csv file (for FlowJo users: select scale values).
# Set working directory:
setwd("path_to_your_folder_with_csv_files")
# Find file names of .csv files in the current working directory:
filenames <- list.files(pattern = ".csv")
# Verify:
filenames

## Defining a function to read a flow cytrometry file in csv format:
# Each row is a cell, each column is a parameter. In our experience, the flow cytometers sometimes output duplicate entries (listing the same cell twice), we remove these and report.
# Please check how your csv file is separated and adjust the sep argument in the function if necessary. In this example we import a semicolon separated file.
read.flow_csv <- function(pathIN){
  raw <- read.csv(pathIN, sep=";", header=TRUE, stringsAsFactors=FALSE)
  IND <- which(duplicated(raw))
  # Check for duplicates and report if found:
  if(any(duplicated(raw))){
    cat(paste0("=== Duplicate entries removed in [",pathIN,"]: ",length(IND)," ===\n"))
    print(head(raw[IND,]))
    cat("----\n")
  }
  return(unique(raw))
}

# Read all:
dfs <- sapply(filenames,read.flow_csv,simplify=FALSE)

##############################
#REWRITE TO FLOWFRAME/FLOWSET# 
##############################

## Defining a function to rewrite a csv into a flowframe:
csv_2_ff <- function(dat){
  # Compute required metadata - column names with description - ranges, min, and max settings
  meta <- data.frame(name=dimnames(dat)[[2]],
                     desc=paste(dimnames(dat)[[2]]),
                     range =(apply(apply(dat,2,range),2,diff)),
                     minRange = apply(dat,2,min),
                     maxRange = apply(dat,2,max))
  # Create flowframe
  flowframef <- new("flowFrame",exprs=as.matrix(dat),parameters=AnnotatedDataFrame(meta))
  return(flowframef)
}

# rewrite to flowframe
dfs_ff = sapply(dfs,function(x) csv_2_ff(x),simplify=FALSE)

# rewrite to flowset
dfs_fs <- as(dfs_ff,"flowSet")


###############################
####### TRANSFORMATION ########
###############################

## Each parameter of interest needs to be arcsinh transformed with an individual cofactor. The cofactor can be deduced from the size of the linear region around zero on a biexponential scale, as plotted in a histogram (in conventional gating software).
# Choose manual transformation or automated transformation (we prefer manual)
# Define parameters and cofactors for transformations:
manualcofactors <- c(CD95=524,CD8=262,CD27=263,CCR7=1787,CD45RA=524,CD3=678,CD49b=398,CD4=915,CD69=830,CD103=504)
automatedcofactors <- estParamFlowVS(dfs_fs, names(manualcofactors)) #this may take a while.

#manual
dfs_fs_t_manual <- transFlowVS(dfs_fs, channels=names(manualcofactors), cofactor=manualcofactors)

#auto
dfs_fs_t_auto <- transFlowVS(dfs_fs, channels=names(manualcofactors), cofactor=automatedcofactors)

#check whether your are satisfied with the transformed parameters
flowViz.par.set(theme =  trellis.par.get(), reset = TRUE)
densityplot(~CD95+CD8+CD27+CCR7+CD45RA+CD3+CD49b+CD4+CD69+CD103, dfs_fs_t_manual, main="manual")
densityplot(~CD95+CD8+CD27+CCR7+CD45RA+CD3+CD49b+CD4+CD69+CD103, dfs_fs_t_auto, main="auto")

##############################
######## NORMALIZATION #######
##############################

## To correct for technical inter-sample variation we apply normalization by fdaNorm (which automatically detects the number of peaks)
# We continue with the manual transformed dataset
# Select the markers which require normalization (based on the densityplots you generated above). Be aware that you don't remove biological variation!
dfs_fs_t_manual_normfda <- warpSet(dfs_fs_t_manual, stains=c('CD8','CD27','CCR7','CD3','CD49b','CD4'))

#check whether you are satisfied
densityplot(~CD8+CD27+CCR7+CD3+CD49b+CD4, dfs_fs_t_manual_normfda, main="fdaNorm")

##############################
####### EXPORT TO FCS ########
##############################

## The flowset (dfs_fs_t_manual_normfda) can be exported to individual fcs files
# Create an 'output' folder
dir.create("Output", showWarnings = FALSE)
setwd("Output")

#Save flowframes wihtin flowset as fcs files using the flowCore package
write.flowSet(dfs_fs_t_manual_normfda, outdir='Output', filename = paste0(gsub(".csv", ".fcs", sampleNames(dfs_fs_t_manual_normfda))))
