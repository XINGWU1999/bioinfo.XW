#' @name GO_plot
#' @title return a GO plot or GO object
#'
#' @details return a GO plot or GO object
#'
#' @param Gene gene to plot
#' @param GO_result a dataframe of GO result
#' @param GO a returned ego from return_ego = T
#' @param input_type type of gene input
#' @param showCategory number of categories to plotted
#' @param ignore_gene defult as F, wether to ignore some GO item related to the character described in ignore_list
#' @param ignore_list default as "RNA", "catabolic", "metabolic", "viral", "translation initiation", "protein target"
#' @param title title to add on GO plot
#' @param label_format default as 30
#' @param font.size default as 12
#' @param return_ego default as F
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
GO_plot <- function(Gene = NULL, GO_result = NULL, GO = NULL, input_type = "SYMBOL", showCategory = 20, title = "",
                    label_format = 30, font.size = 12, return_ego = F, ignore_gene = F,
                    ignore_list = c("RNA", "catabolic", "metabolic", "viral", "translation initiation", "protein target"),
                    Want_gene = F, Want_list = NULL,...){
  GO_S4 = GO
  if(is.null(Gene) & is.null(GO_S4) & !is.null(GO_result)){
    #GO_result present
    if(ignore_gene)  GO_result <- dplyr::filter(GO_result,
                                         !stri_detect_regex(GO_result$Description, stri_join(ignore_list, collapse = "|")))
    GO_S4 = new("enrichResult", result = GO_result)
    if(return_ego) return(GO_S4)
    return(enrichplot::dotplot(GO_S4, showCategory = showCategory, title = title, label_format = label_format, font.size = font.size, x = "Count", ...))
  }


  if(is.null(Gene) & is.null(GO_result) & !is.null(GO_S4)){
    # GO_S4 (EGO result) present
    if(ignore_gene) GO_S4 = new("enrichResult", result = dplyr::filter(GO_S4@result,
                                                                !stri_detect_regex(GO_S4@result$Description, stri_join(ignore_list, collapse = "|"))) )
    if(Want_gene) GO_S4 = new("enrichResult", result = dplyr::filter(GO_S4@result,
                                                              stri_detect_regex(GO_S4@result$Description, stri_join(Want_list, collapse = "|"))) )
    if(return_ego) return(GO_S4)
    return(enrichplot::dotplot(GO_S4, showCategory = showCategory, title = title, label_format = label_format, font.size = font.size, x = "Count", ...))
  }


  if(input_type != "ENTREZID" && input_type != "entrezid") Gene <- convert(ENSG_list = Gene, from_type = input_type, to_type = "ENTREZID")
  # gene present
  library(clusterProfiler)
  library(org.Hs.eg.db)
  ego <- enrichGO(gene = Gene,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE)
  ego_result = dplyr::filter(ego@result, pvalue < 0.05)
  if(ignore_gene) ego = new("enrichResult", result = dplyr::filter(ego_result,
                                                            !stri_detect_regex(ego_result$Description, stri_join(ignore_list, collapse = "|"))) )
  if(Want_gene) ego = new("enrichResult", result = dplyr::filter(ego_result,
                                                            stri_detect_regex(ego_result$Description, stri_join(Want_list, collapse = "|"))) )

  if(return_ego) return(ego)
  return(enrichplot::dotplot(ego, showCategory = showCategory, title = title, label_format = label_format, font.size = font.size, x = "Count", ...))
  detach(clusterProfiler)
}
