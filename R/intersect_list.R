#' @name intersect_list
#' @title return
#'
#' @details
#'
#' @param target_list the list of element to be intersected
#' @export
#'

intersect_list <- function(target_list){
  if(class(target_list) != "list") return(target_list)
  if(length(target_list) == 1) res <- target_list[[1]]
  if(length(target_list) > 1) res <- intersect(target_list[[1]], target_list[[2]])
  if(length(target_list) > 2) for(i in 3:length(target_list)) res <- intersect(res, target_list[[i]])
  return (res)
}
