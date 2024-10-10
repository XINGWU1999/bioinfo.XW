#' @name Dim_red_Cluster
#' @title Dim reduction and find cluster of Seurat object
#'
#' @details to perform dim reduction and cluster manually
#'
#' @param work_object Seurat object
#' @param npcs default as 20
#' @param ... other parameter of FindVariableFeatures
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
Dim_red_Cluster <- function(work_object, npcs = 20, ...){
  cat("Normalize & Scaling")
  work_object_scale <- work_object %>%
    NormalizeData(verbose = F) %>%
    FindVariableFeatures(verbose = F, ...) %>%
    ScaleData(verbose = F)

  cat("Dim reduction")
  work_object_DR <- work_object_scale %>%
    RunPCA(npcs = npcs, features = VariableFeatures(.), verbose = F) %>%
    RunUMAP(dims = 1:npcs, verbose = F)

  #Cluster
  cat("Clustering")
  work_object_Clst <- work_object_DR %>%
    FindNeighbors(dims = 1:npcs)

  for(resi in seq(0,2,0.1)){
    work_object_Clst %<>%
      FindClusters(res = resi, verbose = F)
  }

  return(work_object_Clst)
}
