# NaraR.r
setwd("Documents/NaraR/")
library(maptools)
library(rgdal)
library(rgeos)
library(PBSmapping, quietly=T)


# country list => i.e. for focal area
focal_list = c("Taiwan", "China", "Korea, Republic of", "Japan", "Philippines", "Viet Nam",
               "Lao People's Democratic Republic", "Cambodia","Thailand","Burma","Malaysia")

### 1. INPUT polygon shapes
# TODO: write a function to check the projection of input shp
species <- readShapePoly("data/Fairy-Pitta.shp",
                         proj4string=CRS("+proj=longlat +datum=WGS84"))
## 1.0 Global data input
# world map from
# http://thematicmapping.org/downloads/world_borders.php
ref_area <- readShapePoly("data/TM_WORLD_BORDERS-0.3.shp", 
              proj4string=CRS("+proj=longlat +datum=WGS84"))
#ref_area.attr <- ShowSpatialPointsDataFrame(ref_area)

# biogeographical assessment: biomes
ref_bio <- readShapePoly("data/Asia_biomes_diss.shp", 
              proj4string=CRS("+proj=longlat +datum=WGS84"))
#ref_bio.attr <- ShowSpatialPointsDataFrame(ref_bio)

## 1.1 PROJECTION 
# reproject to cyclindrial equation area EPSG 3410
proj.string <- CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +a=6371228 +b=6371228 +units=km +no_defs")
ref_area.prj <- rgdal::spTransform(ref_area, CRS=proj.string)
ref_bio.prj <- rgdal::spTransform(ref_bio, CRS=proj.string)
species.prj <- rgdal::spTransform(species, CRS=proj.string)



## [FUNCTIONS]
## [FUNCTION] Calculate area for given projected polygon
calc_area <- function(SpatialPolygon.prj) {
  require(PBSmapping)
  SpatialPolygon.ps <- combinePolys(SpatialPolygons2PolySet(SpatialPolygon.prj))
  area <- sum(calcArea(SpatialPolygon.ps, rollup=1)$area) 
  return(area)
}

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

### 2. GLOBAL assessment

# [FUNCTION::Gloabal_assessment] Find biogeographical region numbers 
# within distribution range
find_over_num <- function(target_layer, background_layer) {
  require(sp)
  intersection.lst <- sp::over(target_layer, background_layer, returnList=T)
  return(length(unique(unlist(intersection.lst))))
}

# [FUNCTION::Global_assessment] Asssessment of global distribution
# biome is the number of distributed biomes for target species (integer)
# The thres_lower indicates the lower cut-off value (integer); thres_upper
# means the upper cut-off value (integer). Default lower cut-off value
# is 1; upper cut-off value is 3
glob_eval <- function(biogr, thres_lower=1, thres_upper=3) {
  if ( biogr-round(biogr) != 0 | thres_lower-round(thres_lower) != 0 | thres_upper-round(thres_upper) != 0) { 
    print("Please input non-negative integer") 
  } else {
    if ( biogr <= thres_lower ) {
      sp_distr <- c(0, "Local")
    } else if ( biogr > thres_lower & biogr <= thres_upper  ) {
      sp_distr <- c(1, "Regional")
    } else if ( biogr > thres_upper  ) {
      sp_distr <- c(2, "Global")
    }
    return(sp_distr)
  }
}

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


NaraR <- function(species.prj, ref_area.prj, focal_list, Area_sp_in_ref, Area_ref, thres=2, gthl=1, gthu=3) {
  species.ps <- combinePolys(SpatialPolygons2PolySet(species.prj))
  ass.glob <- glob_eval(find_over_num(species.prj, ref_bio.prj), gthl, gthu)
  ass.reg <- reg_eval(species.ps, ref_area.prj, focal_list, Area_sp_in_ref, Area_ref, thres=2)
  
  glob_a_reg <- as.numeric(ass.glob[1])+ as.numeric(ass.reg[2,])
  glob_a_reg[glob_a_reg < 0]<- "NA"
  glob_a_reg[glob_a_reg == 3] <- "Basic"
  glob_a_reg[glob_a_reg == 2] <- "Medium"
  glob_a_reg[glob_a_reg == 1] <- "High"
  glob_a_reg[glob_a_reg == 0] <- "Very High"
  returnNaraR <- rbind(ass.reg[1,],glob_a_reg)
  return(returnNaraR)
}

# calculate Area of Reference first
ref_area.ps <- combinePolys(SpatialPolygons2PolySet(ref_area.prj))
species.ps <- combinePolys(SpatialPolygons2PolySet(species.prj))
Area_ref <- calc_area(ref_area.prj)
Area_sp_in_ref <- sum(calcArea(combinePolys(joinPolys(species.ps, ref_area.ps, "INT", maxVert=1.0e+06)), rollup=1)$area)

NaraR(species.prj, ref_area.prj, focal_list, Area_sp_in_ref, Area_ref, thres=2, gthl=1, gthu=3)





