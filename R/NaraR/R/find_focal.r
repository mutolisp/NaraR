## [FUNCTION] Find countries in focal list and store PolySet in a vector
find_focal <- function(ref_area.prj, focal_list) {
  require(PBSmapping)
  focal_list.ps <- vector(mode="list", length=length(focal_list)) 
  for ( i in 1:length(focal_list) ) {
    # check if the focal_list names match the fieldname in ref_area.prj 
    if ( focal_list[i] %in% ref_area.prj$NAME ) {
      focal_list.ps[[i]] <- combinePolys(SpatialPolygons2PolySet(ref_area.prj[ref_area.prj$NAME==focal_list[i],]))
    } else { 
      print(paste("Error!", focal_list[i],"does not match in the reference area", sep=" ")) 
    }
  }
  return(focal_list.ps)
}
