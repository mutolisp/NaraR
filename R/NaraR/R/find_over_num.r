# [FUNCTION::Gloabal_assessment] Find biogeographical region numbers 
# within distribution range
find_over_num <- function(target_layer, background_layer) {
  require(sp)
  intersection.lst <- sp::over(target_layer, background_layer, returnList=T)
  return(length(unique(unlist(intersection.lst))))
}
