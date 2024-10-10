#' @name convert
#' @title To convert among gene nomenclature among ENSG, symbol, entrezid .etc
#'
#' @details return either a list of input genes or a corresponding data frame
#'
#' @param ENSG_list input list
#' @param from_type type of input list
#' @param type type of output list
#' @param return_list whether to return a list of gene (or a data frame insteat)
#' @export
#'
#'

library(Seurat)
library(stringi)
library(magrittr)
library(ggplot2)
library(tidyverse)
library(Matrix)

convert <- function(ENSG_list, from_type = "GENEID",  to_type = "SYMBOL", return_list = T, ENSG_provided = NULL, Symbol_provided = NULL){
  mutate_list = function(xlist, x_cond, y_cond){
    res <- c()
    for(i in 1: length(xlist)) res = c(res, y_cond[which(x_cond == xlist[i])])
    return(res)
  }
  library(EnsDb.Hsapiens.v86)
  corr_X = c("ENSG",   "ensg",  "geneid", "Geneid", "GENEID", "symbol","Symbol","name", "SYMBOL", "ENTREZID", "entrezid")
  corr_Y = c("GENEID", "GENEID", "GENEID", "GENEID","GENEID", "SYMBOL","SYMBOL","SYMBOL", "SYMBOL", "ENTREZID", "ENTREZID")
  order = mutate_list(c(from_type, to_type), corr_X, corr_Y)
  if(!is.null(ENSG_provided) & ! is.null(Symbol_provided) & all(order %in% c("SYMBOL", "GENEID"))){
    if(order[[1]] == "GENEID") gene_list = data.frame(SYMBOL = mutate_list(ENSG_list, x_cond = Gene_info$gene_id, y_cond = Gene_info$gene_name),
                                                      GENEID = ENSG_list)
    if(order[[1]] == "SYMBOL") gene_list = data.frame(GENEID = mutate_list(ENSG_list, x_cond = Gene_info$gene_name, y_cond = Gene_info$gene_id),
                                                      SYMBOL = ENSG_list)
  }else gene_list <- ensembldb::select(EnsDb.Hsapiens.v86, keys = ENSG_list, keytype = order[[1]], columns = c(order[[2]], "GENEID"))

  if(return_list){
    pb <- txtProgressBar(style=3)
    res <- c()
    for(i in 1: length(ENSG_list)){
      each = ENSG_list[i]
      setTxtProgressBar(pb, i/length(ENSG_list))
      res_each = ifelse(each %in% gene_list[,order[[1]]],
                        gene_list[which(gene_list[,order[[1]]] == each), order[[2]]],
                        NA)
      res = c(res, res_each)
    }
  }else res = gene_list
  warning(stri_join(round(dim(gene_list)[1]/length(ENSG_list)*100,2) , "% of genes are successfully converted"))
  detach("package:EnsDb.Hsapiens.v86", unload = TRUE)
  detach("package:ensembldb", unload = TRUE)
  return(res)
}
#这是一个很好用的函数，ENSG_list就是一个ENSG开头的geneid向量即可，默认的转换type为“SYMBOL”
#可选的转换方向包括"GENEID", "GENENAME", "TXID", "TXBIOTYPE", "TXSEQSTART","TXSEQEND", "SEQNAME", "SEQSTRAND"
#一次只能选一个
