#' @name header
#' @title glance at the data
#'
#' @details return 8X8 if dimension < 8X10
#'
#' @param df the 2d-matrix or dataframe
#' @param x dimension 1 to be returned
#' @param y dimension 2 to be returned
#' @export
#'
#'
#'
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
