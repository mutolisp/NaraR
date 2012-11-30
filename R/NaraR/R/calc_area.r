## [FUNCTION] Calculate area for given projected polygon
calc_area <- function(SpatialPolygon.prj) {
  require(PBSmapping)
  SpatialPolygon.ps <- combinePolys(SpatialPolygons2PolySet(SpatialPolygon.prj))
  area <- sum(calcArea(SpatialPolygon.ps, rollup=1)$area) 
  return(area)
}

