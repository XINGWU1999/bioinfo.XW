#' @name SearchPM
#' @title search PM result
#'
#' @details return a data.frame, better to check the page numbers (using return URL = T) before searching
#'
#' @param keyword search the keyword
#' @param start_pages num of page start to search
#' @param num_pages end of page to search
#' @param years a contiuous year, like c(2022:2024)
#' @param dir_save dir to save
#' @param returnres return the result apart from saving the results
#' @param returnURL return the URL to check the page nums
#'
#'
#'
# 导入必要的库
library(rvest)
library(stringi)
library(tidyverse)
#总体函数
SearchPM = function(keyword = NULL, start_pages = 1, num_pages = 1, years = c(2020:2024), dir_save = "./", returnres = F, returnURL = F){
  # 设置关键词和页数
  # keyword <- "continuous glucose monitor"
  # num_pages <- 5  # 假设我们要提取前5页的信息

  if(returnURL){
    URL = paste0("https://pubmed.ncbi.nlm.nih.gov/?term=", URLencode(keyword), "&page=", num_pages, "&format=abstract&size=10")
    if(!is.null(years)) URL = paste0(URL, "&filter=years.", min(years), "-", max(years))
    return(URL)
  }
  #将空结果转化为NA的函数
  convert_empty_to_na <- function(obj) {
    if(length(obj) == 0) {
      return(NA)
    } else {
      return(obj)
    }
  }
  # 定义函数以获取目标信息
  get_pubmed_info <- function(keyword, page_num, year = years) {
    # 构建搜索URL
    search_url <- paste0("https://pubmed.ncbi.nlm.nih.gov/?term=", URLencode(keyword), "&page=", page_num, "&format=abstract&size=10")
    if(!is.null(year)) search_url = paste0(search_url, "&filter=years.", min(year), "-", max(year))

    # 抓取网页内容
    page <- read_html(search_url)
    # 储存剩余的搜索结果数目
    All_res_num = page %>% html_nodes(xpath = stri_join('/html/body/main/div[9]/div[2]/div[2]/div[1]/div[1]/h3/span')) %>%
      html_text() %>% stri_replace_all_fixed(",", "") %>% as.numeric()
    First_num =  page %>% html_nodes(xpath = stri_join("/html/body/main/div[9]/div[2]/section[1]/div[2]/div/div[2]/div/label/span")) %>%
      html_text() %>% stri_replace_all_fixed(",", "") %>% as.numeric()

    # 提取信息
    info <- data.frame()
    for (i in 1:min(10, All_res_num-First_num + 1)) { # 提取每页的前五条信息
      # 提取信息
      pmid <- page %>% html_nodes(xpath = stri_join('/html/body/main/div[9]/div[2]/section[1]/div[2]/div/div[', i, ']/article/header/div[1]/ul/li[1]/span/strong')) %>% html_text() %>% trimws() %>% convert_empty_to_na
      abstract <- page %>% html_nodes(xpath = paste0('/html/body/main/div[9]/div[2]/section[1]/div[2]/div/div[', i, ']/article/div[1]')) %>% html_text() %>% trimws() %>% convert_empty_to_na
      title <- page %>% html_nodes(xpath = paste0('/html/body/main/div[9]/div[2]/section[1]/div[2]/div/div[', i, ']/article/header/div[1]/h1/a')) %>% html_text() %>% trimws() %>% convert_empty_to_na
      journal <- page %>% html_nodes(xpath = paste0('/html/body/main/div[9]/div[2]/section[1]/div[2]/div/div[', i, ']/article/header/div[1]/div[1]/div[2]/div/button')) %>% html_text() %>% trimws() %>% convert_empty_to_na


      # 组合信息到data.frame
      info <- rbind(info, data.frame(PMID = pmid, Abstract = abstract, Title = title, Journal = journal))
    }

    return(info)
  }
  # 定义循环以提取每一页的信息

  pb <- txtProgressBar(style=3) #设置一个进度条 #循环前

  #循环提取每一页的函数
  scrape_pubmed_pages <- function(keyword, num_pages) {
    all_info <- data.frame()
    for (i in start_pages:num_pages) {
      setTxtProgressBar(pb, i/num_pages) #循环内
      info <- get_pubmed_info(keyword, i)
      all_info <- rbind(all_info, info)
    }
    return(all_info)
  }

  all_info <- scrape_pubmed_pages(keyword, num_pages) %>% as.data.frame() %>% filter(!is.na(Title))
  save(all_info, file = paste0(dir_save, "/", keyword, "_", min(years), "-", max(years), "res.RData"))
  write.csv(all_info, paste0(dir_save, "/", keyword, "_", min(years), "-", max(years), "res.csv"))
  if(returnres) return(all_info)
}
#
# Res = SearchPM(keyword = "continuous glucose monitor", num_pages = 5, years = c(2020:2024), dir_save = "./test", returnURL = F)
