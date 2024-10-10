#' @name daily
#' @title generate daily file fastly and easily
#'
#' @details just generate 2 files (data and scripe)
#'
#' @param date date
#' @param name name or title of today's work
#' @param project which project is the work is under, currently available (Module_correlation, Remove_batch, Dist_algorithm, data_generation)
#' @export
#'
#'
library(stringi)
daily <- function(date, name, project){
  path = switch (project,
    "Module_correlation" = "~/Desktop/Lab/Na_Lab/Project/CMO_seq_analysis/Cross_analysis_preprocess/Module_correlation/",
    "Remove_batch" = ".~/Desktop/Lab/Na_Lab/Project/CMO_seq_analysis/Cross_analysis_preprocess/Remove_batch_correlation/",
    "Dist_algorithm" = "~/Desktop/Lab/Na_Lab/Project/CMO_seq_analysis/Cross_analysis_preprocess/Dist_Algorithm/",
    "data_generation" = "~/Desktop/Lab/Na_Lab/Project/MSC-tool/Data_generation/"
  )

  dir.create(stri_join(path, date,"_", name))
  dir.create(stri_join(path, date,"_", name, "/script"))
  dir.create(stri_join(path, date,"_", name, "/data"))

}
