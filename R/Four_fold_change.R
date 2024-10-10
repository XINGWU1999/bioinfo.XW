#' @name Four_fold_change
#' @title add Sensitivity, NPV, accuracy, etc. following the four-fold table
#'
#' @details add Sensitivity, NPV, accuracy, etc. following the four-fold table
#'
#' @param Four_fold a dataframe, colnames as(c("Predict", "True")), Predict has c("Negative", "Positive"), True has c("Ctrl", "Disease")
#' @export
#'
#'
#'
library(magrittr)
library(tidyverse)
Four_fold_change = function(Four_fold){
  Four_fold = Four_fold %>%
    cbind(apply(., 1, sum)) %>%
    rbind(apply(., 2, sum)) %>%
    set_colnames(c(colnames(.)[1:2], "Sum")) %>%
    set_rownames(c(rownames(.)[1:2], "Sum"))
  New_Four_fold = rbind(Four_fold,
                        data.frame(Ctrl = "", Disease = "", Sum = "", row.names = ""),
                        data.frame(Ctrl = Four_fold[2,2]/Four_fold[3,2], Disease = "", Sum = "", row.names = "Sensitivity"),
                        data.frame(Ctrl = Four_fold[1,1]/Four_fold[3,1] , Disease = "", Sum = "", row.names = "Specificity"),
                        data.frame(Ctrl = Four_fold[2,2]/Four_fold[2,3] , Disease = "", Sum = "", row.names = "PPV"),
                        data.frame(Ctrl = Four_fold[1,1]/Four_fold[1,3], Disease = "", Sum = "", row.names = "NPV"),
                        data.frame(Ctrl = ((Four_fold[1,1] + Four_fold[2,2])/Four_fold[3,3]), Disease = "", Sum = "", row.names = "Accuracy")
  )
  return(New_Four_fold)
}
