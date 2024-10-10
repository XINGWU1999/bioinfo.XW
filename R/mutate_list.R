#' @name mutate_list
#' @title convert list x into y through the map relationship x_cond --> y_cond
#'
#' @details return a list as long as the length of x_cond = y_cond, if a element is present in x_list but absent in x_cond, return NA in y_list
#'
#'
#' @param xlist xlist
#' @param x_cond x_cond
#' @param y_cond y_cond
#' @param force_1st if there are some replication in x_cond, whether to use the 1st element occured in x_cond to match
#' @param use_origin if one element in xlist is not present in x_cond, then use the origin element in xlist(T) or use NA(F), default as F
#' @export
#'
#'

library(stringi)
library(magrittr)
library(ggplot2)
library(tidyverse)
library(Matrix)
mutate_list = function(xlist, x_cond, y_cond, force_1st = F, use_origin = F, verbose = T){
  if(length(x_cond) != length(y_cond)){
   return(0)
    warning("Length of X and Y are not matched")
  }#不一一对应直接返回0

  if(!all(unique(xlist) %in% x_cond)){
    if(!use_origin) warning("There will be some NA in y_list")
  }

  #x_list 内有不在x_cond 内的元素，使用NA代替
  if(force_1st){ #如果确定强制转换（使用第一个）
    if(length(x_cond) != length(unique(x_cond))) warning("There are replications in x_cond, will automatically choose the 1st x_cond")
    #x_cond 内有重复元素，说明存在多对1的映射关系，可选是否直接选择第一个
    res <- c()
    if(verbose) pb <- txtProgressBar(style=3)
    for(i in 1: length(xlist)) {#循环加入元素
      res = c(res, ifelse(xlist[i] %in% x_cond, #判断该元素是否在x_cond内
                          ifelse(length(which(x_cond == xlist[i])) > 1, #判断该元素在x_cond内是否出现多次
                                 y_cond[which(x_cond == xlist[i])[1]],  #如果出现多次则返回y_cond内对应的第一个元素
                                 y_cond[which(x_cond == xlist[i])]),   #如果仅一次则直接返回
                          ifelse(use_origin, xlist[i], NA) #如果不在直接返回NA
                          ))
      if(verbose) setTxtProgressBar(pb, i/length(xlist))
    }
  }else{#如果不强制转换，则当有重复元素时直接返回0
    if(length(x_cond) != length(unique(x_cond))) {
      warning("replicated data in x_cond")
      return(0)
      }

    #x_cond 内有重复元素，说明存在多对1的映射关系，可选是否直接选择第一个
    if(verbose) pb <- txtProgressBar(style=3)
    res <- c()
    for(i in 1: length(xlist)) {
      res = c(res, ifelse(xlist[i] %in% x_cond,
                          y_cond[which(x_cond == xlist[i])],
                          ifelse(use_origin, xlist[i], NA) #如果不在直接返回NA
                          ))
      if(verbose) setTxtProgressBar(pb, i/length(xlist))
    }
  }
  return(res)
}
