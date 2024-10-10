#' @name DEG_generation
#' @title generate differential expressed gene
#'
#' @details
#'
#'
#' @param work_object
#' @param optmz_res resolution chosen for DEG genetion
#' @param p_thr p value threshold
#' @param logFC_thr logFC threshold
#' @param test_use default "wilcox" details from FindMarkers
#' @param only.pos only return the positions
#' @export
#'
#'

library(Seurat)
library(stringi)
library(magrittr)
library(ggplot2)
library(tidyverse)
library(Matrix)
DEG_generation <- function(work_object, optmz_res = NA,
                           p_thr = 0.05, logFC_thr = 0.25, test_use = "wilcox",
                           only.pos = T, ...){
  p_value_score <- function(p_val){
    res = c()
    for(each in p_val){
      res = c(res, ifelse(each < 0.05,
                          -log2(each),
                          ifelse(each >= 0.1,
                                 0,
                                 -86.4385*each + 8.643)))
    }
    return(res)
  }
  if(is.na(optmz_res)) return(NA)
  Idents(work_object) <- stri_join("RNA_snn_res.",optmz_res)
  Marker_list = list()
  for(i in 1:length(unique(Idents(work_object)))){
    each_cluster = sort(unique(Idents(work_object)))[i]
    Marker_list[[i]] <- FindMarkers(work_object, ident.1 = each_cluster,
                                    return.thresh = p_thr, logfc.threshold = logFC_thr,only.pos = only.pos, ...) %>%
      mutate(p_val_score = p_value_score(p_val))
  }
  return(Marker_list)
}
