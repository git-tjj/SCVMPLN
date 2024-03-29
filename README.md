# SCVMPLN: Variational inference for Mixture Poisson Log-Normal graphial model in Single-cell gene regulatory network analysis with mixed cell populations
> SCVMPLN is a R package for cell-type-specifc gene regulatory inference for Single-cell RNA-seq data with mixed cell populations using a block-wise descent algorithm based on the variational inference .

## Installation
**SCVMPLN** is available on [Github](https://github.com/git-tjj/SCVMPLN).

### R Package installation
- For the development version, use the github install
```{r package github, eval = FALSE}
remotes::install_github("https://github.com/git-tjj/SCVMPLN")
```

## Usage and main functions
The package comes with a simulation example single-cell RNA-seq data to present the functionality of main function.

### Preparsion
```{r load SCVMPLN, eval = FALSE}

##1. Package loading
##-------------------------------------------------------------
library(SCVMPLN)
library(Seurat)
##-------------------------------------------------------------

##2. load reference dataset
##-------------------------------------------------------------
data("example_data")
##-------------------------------------------------------------

```

### Initialization
```{r, warning = FALSE}
##-------------------------------------------------------------
gene_GRN <- rownames(example_data)
VMPLN_list<-VMPLN_init(expression_profile = example_data,
                     celltypes_num = 3,
                     celltypes_ref = NULL,
                     ls_est = "TSS",
                     gene_GRN = gene_GRN,
                     HVG_model_num = 0,
                     zero_GRN = NULL,
                     preprocess_Control = list(HVG_num = length(gene_GRN),npc = 50,
                     run_umap = TRUE,label_umap = NULL,
                     cluster_method = "Kmeans",resolution = 0.8),
                     core_num = 1)
##-------------------------------------------------------------
```
### (Optional) Exact the path of hyper-parameter lambda
```{r, warning = FALSE}
##-------------------------------------------------------------
gene_GRN <- rownames(example_data)
VMPLN_list<-VMPLN_init(expression_profile = example_data,
                     celltypes_num = 3,
                     celltypes_ref = NULL,
                     ls_est = "TSS",
                     gene_GRN = gene_GRN,
                     HVG_model_num = 0,
                     zero_GRN = NULL,
                     preprocess_Control = list(HVG_num = length(gene_GRN),npc = 50,
                     run_umap = TRUE,label_umap = NULL,
                     cluster_method = "Kmeans",resolution = 0.8),
                     core_num = 1)
##-------------------------------------------------------------
```
###  Run the main function for optimization of model
```{r, warning = FALSE}
##-------------------------------------------------------------
VMPLN_list_alongpath<-list()
for(l in 1:length(lambda_path)){
cat(paste(paste(rep("##",40),sep = "",collapse = ""),"\n",(paste(rep("##",40),sep = "",collapse = "")),"\n",sep = "",collapse = ""))
##
VMPLN_list_alongpath[[l]]<-VMPLN_main(VMPLN_list = VMPLN_list,lambda_use = lambda_path[l],
                                      Theta_Control = list(penalize_diagonal = FALSE),
                                      U_fix = FALSE,
                                      verbose = TRUE,
                                      core_num = 1)
}
names(VMPLN_list_alongpath)<-paste("lambda: ",lambda_path,sep = "")
##-------------------------------------------------------------
```

## Reference

Please cite our work using the following references:

Tang, J., Wang, C., Xiao, F., & Xi, R. (2022). Network Analysis of Count Data from Mixed Popula+ons. arXiv preprint arXiv:2212.03665.
