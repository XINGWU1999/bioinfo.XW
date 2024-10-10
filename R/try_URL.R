#' @name try_URL
#' @title try to open a URL link
#'
#' @details return read_html(URL)
#'
#' @param URL the link
#' @param upperlimit max time to try
#' @param content expressed on the screen
#' @export
#'
#'
#'
library(xml2)
try_URL <- function(URL, upperlimit = 3, content = ""){
  web = NA
  i = 1
  while (is.na(web) & i <= upperlimit){
    print(stri_join("trying to connect to NCBI for ", content, "(time:", i, ")"))
    web <- tryCatch({
      read_html(URL)
    },#这个逗号千万别删了。。。不然报错
    error = function(e){
      cat("ERROR :",conditionMessage(e),"\n")
      return(NA)})
    i <- i+1
  }

  if(i >= upperlimit & is.na(web)){
    warning("bad connection")
    return(NA)
  }

  return(web)
}
