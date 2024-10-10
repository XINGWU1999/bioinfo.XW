#' @name QC_subset
#' @title return a subset of Surat object
#'
#' @details return 8X8 if dimension < 8X10
#'
#' @param Count_low
#' @param Count_high
#' @param Feature_low
#' @param Feature_high
#' @param percent.mt_high default = 100, no filter on percent of mt gene
#' @param Assign whether to assign the name of cell or gene, if T, then cell_assign_target and gene_assign_target should be assigned
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
QC_subset <- function(work_object,
                      Count_low = NA, Count_high = NA,
                      Feature_low = NA, Feature_high = NA,percent.mt_high = 100,
                      Assign = F, cell_assign_target = NA, gene_assign_target = NA){
  work_object[["percent.mt"]] <- PercentageFeatureSet(work_object, pattern = "^MT-")
  if(Assign & (! any(is.na(cell_assign_target))) & (!any(is.na(gene_assign_target)))) return(work_object %>%
                                                                                               subset(subset = cell.name %in% cell_assign_target,
                                                                                                      features = gene_assign_target))
  if(!is.na(sum(Count_low, Count_high, Feature_low, Feature_high, percent.mt_high))) return(work_object %>%
                                                                                              subset(subset = nFeature_RNA > Feature_low & nFeature_RNA < Feature_high &
                                                                                                       nCount_RNA > Count_low & nCount_RNA < Count_high &
                                                                                                       percent.mt < percent.mt_high))
  return(0)
}
