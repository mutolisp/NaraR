# [FUNCTION::Regional_assessment] Assessment of regional distribution
# require calc_area()
reg_eval <- function(species.ps, ref_area.prj, focal_list, Area_sp_in_ref, Area_ref, thres=2) {
  require(PBSmapping)
  Area_focal <- matrix("NA", nrow=1, length(focal_list))
  Area_sp_in_focal <- matrix("NA", nrow=1, length(focal_list))
  for ( i in 1:length(focal_list) ) {
    print(paste("calculating",focal_list[i], sep=" "))
    Area_focal[i] <- calc_area(ref_area.prj[ref_area.prj$NAME==focal_list[i],])
    try_result <- try(SpatialPolygons2PolySet(ref_area.prj[ref_area.prj$NAME==focal_list[i],]), silent=T)
    #####
    if ( class(try_result)[1] == "try-error" ) {
      Area_sp_in_focal[i] <- "NA"
      print("reg_eval(), converting polygon to polyset error")
    } else {
      focal.ps <- combinePolys(SpatialPolygons2PolySet(ref_area.prj[ref_area.prj$NAME==focal_list[i],]))
      # check if the joinPolys operates well
      try_joinPolys <- try(joinPolys(species.ps, focal.ps, "INT", maxVert=1.0e+06), silent=T)
      if ( class(try_joinPolys)[1] == "try-error") {
        Area_sp_in_focal[i] <- "NA"
      } # if the join returns NULL (no intersection), assign 0 to Area_sp_in_focal
      else if ( class(try_joinPolys)[1] == "NULL" ) {
        Area_sp_in_focal[i] <- 0
      } else {
        Area_sp_in_focal[i] <- sum(calcArea(combinePolys(joinPolys(species.ps, focal.ps, "INT", maxVert=1.0e+06)), rollup=1)$area)
      }
    }
    ######
  }
  # Calculate Distribution of observed reference area
  DPob <- Area_sp_in_ref/Area_ref
  # Calculate distribution of expected focal area
  DPexp <- as.numeric(Area_sp_in_focal)/as.numeric(Area_focal)
  reg_ass <- DPob/DPexp
  reg_ass[which(reg_ass <= thres )] <- 1
  reg_ass[which(reg_ass > thres & reg_ass < Inf)] <- 0
  reg_ass[which(reg_ass == Inf)] <- "-9999"
  return_Matrix <- rbind(Area_sp_in_focal, as.numeric(reg_ass))
  colnames(return_Matrix) <- c(focal_list)
  return(return_Matrix)  
}
