#' @name QC_plot
#' @title plot on QC quality
#'
#' @details plot on QC quality
#'
#' @param input can be count matrix or Seurat object
#' @param project_name specify a project_name if input a count matrix
#' @param group use one group to split the qc plot
#' @param organism default to "hs"; can also be set to "ms", so that also mitochondrial genes are filtered as "^mt-"
#' @export
#'
#'
#'

library(Seurat)
library(stringi)
library(magrittr)
library(ggplot2)
library(tidyverse)
library(Matrix)
QC_plot <- function(input, project_name = "test_data", group = NULL, organism = "hs"){
  #Create Seurat
  if(class(input) %in% c("Matrix", "dgCMatrix")){
    work_object <- CreateSeuratObject(counts = input, project = project_name, min.cells = 1, min.features = 50)
  }
  if(class(input) %in% c("Seurat", "SeuratObject")) {
    work_object <- input
  }
  if(organism == "hs") work_object[["percent.mt"]] <- PercentageFeatureSet(work_object, pattern = "^MT-")
  if(organism == "ms") work_object[["percent.mt"]] <- PercentageFeatureSet(work_object, pattern = "^mt-")
  #QC_plot
  P_saturation <- FeatureScatter(work_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = group)

  qtl_mark <- c("99%", "95%", "mid", "5%", "1%")
  nCount_qtl <- round(quantile(work_object@meta.data$nCount_RNA, probs = c(0.99, 0.95, 0.5, 0.05, 0.01)),2)
  nFeature_qtl <- round(quantile(work_object@meta.data$nFeature_RNA, probs = c(0.99, 0.95, 0.5, 0.05, 0.01)),2)
  percent.mt_qtl <- round(quantile(work_object@meta.data$percent.mt, probs = c(0.99, 0.95, 0.5, 0.05, 0.01)), 2)
  P_QC_nCount <- VlnPlot(work_object, features = c("nCount_RNA"), group.by = group) +
    labs(tag = paste(stri_join(qtl_mark, ": ", nCount_qtl), collapse = "\n")) +
    theme(plot.tag.position = "right",
          legend.position = "none",
          plot.tag = element_text(size = 8))
  P_QC_nFeature <- VlnPlot(work_object, features = c("nFeature_RNA"), group.by = group) +
    labs(tag = paste(stri_join(qtl_mark, ": ", nFeature_qtl), collapse = "\n")) +
    theme(plot.tag.position = "right",
          legend.position = "none",
          plot.tag = element_text(size = 8))

  P_QC_percent.mt <- VlnPlot(work_object, features = c("percent.mt"), group.by = group) +
    labs(tag = paste(stri_join(qtl_mark, ": ", percent.mt_qtl), collapse = "\n")) +
    theme(plot.tag.position = "right",
          legend.position = "none",
          plot.tag = element_text(size = 8))

  return(cowplot::plot_grid(P_QC_nCount, P_QC_nFeature, P_QC_percent.mt, P_saturation))
}
