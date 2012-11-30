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
