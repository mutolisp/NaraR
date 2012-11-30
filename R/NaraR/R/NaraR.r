## NaraR
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
