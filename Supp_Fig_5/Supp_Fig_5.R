library(Seurat)
library(SeuratDisk)

# Convert AnnData object to seurat object using SeuratDisk package

#https://github.com/mojaveazure/seurat-disk/issues/109

Convert("./Wimmers/CD14_mono_wimmers_vs_pCoV_12132024.h5ad", dest = "h5seurat", assay = "RNA", overwrite = F, verbose = T)
obj <- LoadH5Seurat("./CD14_mono_cleaned_pCov40_11232021.h5seurat", assays = "RNA")
dim(obj@assays$RNA)

obj <- LoadH5Seurat("./Wimmers/CD14_mono_wimmers_vs_pCoV_12132024.h5seurat",  meta.data = FALSE, misc = FALSE)
print(obj)

#' convertToPseudoBulk function
#'
#' Returns a matrix of pseudobulk expression from a given Seurat object with the selected grouping. 
#' Optionally returns a metadata dataframe fro the samples extracted from the Seurat object. 
#'
#' @param pbmc a seurat object
#' @param groups a character string of variables used for grouping the samples for aggregation
#' @param method method used for aggregating the samples. Can be either one of 'agg'(DEFAULT) or 'avg'. 
#' If 'agg', the expression of cells in the grouping are summed up. If 'avg', the average of the expression of cells in the grouping is used. 
#' @param assays assay from the seurat object to get the expression matrix. 
#' @param slot slot of the assay in the seurat object to use. Can be 'counts'(DEFAULT) for raw data or 'data' for log-normalized expression.
#' @param return.metadata If TRUE(DEFAULT), return the metadata of the samples extracted from the seurat object based on grouping along with the pseudobulk matrix. 
#' If FALSE, only return the pseudobulk matrix.
#' 
#' @return pseudobulk matrix with the provided grouping. Optionally a dataframe of metadata for the samples.
#' @export
#' 
#' 

convertToPseudoBulk <- 
  function (pbmc, 
            groups,
            method = "agg",
            assays = "RNA",
            slot = "counts",
            return.metadata = T){
    
    if(method == "agg"){
      pbmc.agg <- Seurat::AggregateExpression(pbmc, return.seurat=F, group.by=groups, assays=assays, slot = slot)
    }else if(method == "avg"){
      pbmc.agg <- Seurat::AverageExpression(pbmc, return.seurat=F, group.by=groups, assays=assays, slot = slot)
    }else{
      stop("No valid aggregation method suplied. Use either 'agg' for sum or 'avg' for average expression of samples.")
    }
    
    if(return.metadata){
      
      pbmc.agg.metadata <- as.data.frame(pbmc@meta.data[!duplicated(pbmc@meta.data[,groups]),])
      rownames(pbmc.agg.metadata) <- do.call(paste, c(pbmc.agg.metadata[c(groups)], sep="_"))
      
      return(list(Aggregate.Matrix = pbmc.agg[[1]],
                  Aggregate.Metadata = pbmc.agg.metadata))
    }
    
    return(pbmc.agg[[1]])
    
    
  }

# Pseubulk data (Per sample or fig_ids)
obj_agg <- convertToPseudoBulk(seu_obj, groups = c("fig_ids")) #"clinical_groups",

head(obj_agg[['Aggregate.Metadata']])




