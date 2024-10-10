#' @name Cluster_plot
#' @title plot on different resolution
#'
#' @details
#'
#' @param work_object_Clst object with cluster
#' @param quiery_res resolution to be plotted
#' @param label whether to add labels
#' @export
#'
#'

library(Seurat)
library(stringi)
library(magrittr)
library(ggplot2)
library(tidyverse)
library(Matrix)
Cluster_plot <- function(work_object_Clst, quiery_res = c(0.1, 0.5, 1, 1.5), label = T, ...){
  P_Clst = list()
  for(i in 1:length(quiery_res)) P_Clst[[i]] = DimPlot(work_object_Clst, group.by = stri_join("RNA_snn_res.",quiery_res[i]), label = label, ...)
  P_Clst <- cowplot::plot_grid(plotlist = P_Clst)
  return(P_Clst)
}
