#' @name Create_object
#' @title create a seurat object from count matrix
#'
#' @details
#'
#' @param project_name name of project
#' @param filter one of auto, manual, no(NO), 0; auto
#' @param ... other parameter of CreatSuratObject, min.cell, min.features
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
fi <- function(Count_data = NA, project_name = "test_data", filter = "auto", ...){
  if(!is.null(dim(Count_data))){
    work_object = switch (filter,
                          "auto" = CreateSeuratObject(counts = Count_data, project = project_name,
                                                      min.cells = 1, min.features = 100),
                          "manual" = CreateSeuratObject(counts = Count_data, project = project_name,
                                                        ...),
                          "no" = CreateSeuratObject(counts = Count_data, project = project_name,
                                                    min.cell = 1, min.features = 1),
                          "No" = CreateSeuratObject(counts = Count_data, project = project_name,
                                                    min.cell = 1, min.features = 1),
                          "0" = CreateSeuratObject(counts = Count_data, project = project_name,
                                                   min.cell = 0, min.features = 0)
    )

    work_object[["percent.mt"]] <- PercentageFeatureSet(work_object, pattern = "^MT-")
    work_object[["cell.name"]] <- rownames(work_object@meta.data)
  }
  return(work_object)
}
