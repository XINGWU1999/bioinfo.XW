#all

library(Seurat)
library(stringi)
library(magrittr)
library(ggplot2)
library(tidyverse)
library(Matrix)
header = function(df, x = 8, y = 8){
  dimension = dim(df)
  if(is.null(dimension)) {
    warning(paste("length of the vector:", length(df)))
    if(length(df) <= 8) x = length(df)
    return(df[1:x])
  }
  if(dimension[1] <= 8) x = dimension[1]
  if(dimension[2] <= 10) y = dimension[2]
  warning(paste("dimension = ",dimension[1], "X", dimension[2]))
  return(df[1:x, 1:y])
}

mutate_list = function(xlist, x_cond, y_cond, force_1st = F){
  if(length(x_cond) != length(y_cond)) return(0)
  #不一一对应直接返回0
  if(!all(unique(xlist) %in% x_cond)) warning("There will be some NA in y_list")
  #x_list 内有不在x_cond 内的元素，使用NA代替
  if(force_1st){ #如果确定强制转换（使用第一个）
    if(length(x_cond) != length(unique(x_cond))) warning("There are replications in x_cond, will automatically choose the 1st x_cond")
    #x_cond 内有重复元素，说明存在多对1的映射关系，可选是否直接选择第一个
    res <- c()
    pb <- txtProgressBar(style=3)
    for(i in 1: length(xlist)) {#循环加入元素
      res = c(res, ifelse(xlist[i] %in% x_cond, #判断该元素是否在x_cond内
                          ifelse(length(which(x_cond == xlist[i])) > 1, #判断该元素在x_cond内是否出现多次
                                 y_cond[which(x_cond == xlist[i])[1]],  #如果出现多次则返回y_cond内对应的第一个元素
                                 y_cond[which(x_cond == xlist[i])]),   #如果仅一次则直接返回
                          NA))                  #如果不在直接返回NA
      setTxtProgressBar(pb, i/length(xlist))
    }
  }else{#如果不强制转换，则当有重复元素时直接返回0
    if(length(x_cond) != length(unique(x_cond))) return(0)
    #x_cond 内有重复元素，说明存在多对1的映射关系，可选是否直接选择第一个
    pb <- txtProgressBar(style=3)
    res <- c()
    for(i in 1: length(xlist)) {
      res = c(res, ifelse(xlist[i] %in% x_cond, y_cond[which(x_cond == xlist[i])], NA))
      setTxtProgressBar(pb, i/length(xlist))
    }
  }
  return(res)
}

convert <- function(ENSG_list, from_type = "GENEID",  type = "SYMBOL", return_list = T){
  mutate_list = function(xlist, x_cond, y_cond){
    res <- c()
    for(i in 1: length(xlist)) res = c(res, y_cond[which(x_cond == xlist[i])])
    return(res)
  }
  library(EnsDb.Hsapiens.v79)
  corr_X = c("ENSG",   "ensg",  "geneid", "Geneid", "GENEID", "symbol","Symbol","name", "SYMBOL", "ENTREZID", "entrezid")
  corr_Y = c("GENEID", "GENEID", "GENEID", "GENEID","GENEID", "SYMBOL","SYMBOL","SYMBOL", "SYMBOL", "ENTREZID", "ENTREZID")
  order = mutate_list(c(from_type, type), corr_X, corr_Y)
  gene_list <- ensembldb::select(EnsDb.Hsapiens.v79, keys = ENSG_list, keytype = order[[1]], columns = c(order[[2]], "GENEID"))
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
  detach("package:EnsDb.Hsapiens.v79", unload = TRUE)
  detach("package:ensembldb", unload = TRUE)
  return(res)
}

GO_plot <- function(Gene = NULL, GO_result = NULL, GO = NULL, input_type = "SYMBOL", showCategory = 20, title = "",
                    label_format = 30, font.size = 12, return_ego = F, ...){
  GO_S4 = GO
  if(is.null(Gene) & is.null(GO_S4) & !is.null(GO_result)){
    GO_S4 = new("enrichResult", result = GO_result)
    if(return_ego) return(GO_S4)
    return(enrichplot::dotplot(GO_S4, showCategory = showCategory, title = title, label_format = label_format, font.size = font.size, x = "Count", ...))
  }
  if(is.null(Gene) & is.null(GO_result) & !is.null(GO_S4)){
    if(return_ego) return(GO_S4)
    return(enrichplot::dotplot(GO_S4, showCategory = showCategory, title = title, label_format = label_format, font.size = font.size, x = "Count", ...))
  }
  if(input_type != "ENTREZID" && input_type != "entrezid") Gene <- convert(ENSG_list = Gene, from_type = input_type, type = "ENTREZID")
  library(clusterProfiler)
  library(org.Hs.eg.db)
  ego <- enrichGO(gene = Gene,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE)
  if(return_ego) return(ego)
  return(enrichplot::dotplot(ego, showCategory = showCategory, title = title, label_format = label_format, font.size = font.size, x = "Count", ...))
  detach(clusterProfiler)
}
plot_grid <- function(...){
  return(cowplot::plot_grid(...))
}

list_combind <- function(list_to_combine){
  for(each in list_to_combine) res <- c(res, each)
  return(res)
}
intersect_list <- function(target_list){
  if(length(target_list) > 1) res <- intersect(target_list[[1]], target_list[[2]])
  if(length(target_list) > 2) for(i in 3:length(target_list)) res <- intersect(res, target_list[[i]])
  return (res)
}

{
  {
    Create_object <- function(Count_data = NA, project_name = "test_data", filter = "auto", ...){
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
    QC_subset <- function(work_object,
                          Count_low = NA, Count_high = NA,
                          Feature_low = NA, Feature_high = NA,percent.mt_high = 100,
                          Assign = F, cell_assign_target = NA, gene_assign_target = NA){
      if(Assign & (! any(is.na(cell_assign_target))) & (!any(is.na(gene_assign_target)))) return(work_object %>%
                                                                                                   subset(subset = cell.name %in% cell_assign_target,
                                                                                                          features = gene_assign_target))
      if(!is.na(sum(Count_low, Count_high, Feature_low, Feature_high, percent.mt_high))) return(work_object %>%
                                                                                                  subset(subset = nFeature_RNA > Feature_low & nFeature_RNA < Feature_high &
                                                                                                           nCount_RNA > Count_low & nCount_RNA < Count_high &
                                                                                                           percent.mt < percent.mt_high))
      return(0)
    }
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
  }

  #Plot functions
  {
    QC_plot <- function(input, project_name = "test_data", group = NULL){
      #Create Seurat
      if(class(input) %in% c("Matrix", "dgCMatrix")){
        work_object <- CreateSeuratObject(counts = input, project = project_name, min.cells = 1, min.features = 50)
        work_object[["percent.mt"]] <- PercentageFeatureSet(work_object, pattern = "^MT-")
      }
      if(class(input) %in% c("Seurat", "SeuratObject")) work_object <- input
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
    Cluster_plot <- function(work_object_Clst, quiery_res = c(0.1, 0.5, 1, 1.5), label = T, ...){
      P_Clst = list()
      for(i in 1:length(quiery_res)) P_Clst[[i]] = DimPlot(work_object_Clst, group.by = stri_join("RNA_snn_res.",quiery_res[i]), label = label, ...)
      P_Clst <- cowplot::plot_grid(plotlist = P_Clst)
      return(P_Clst)
    }
  }

}

