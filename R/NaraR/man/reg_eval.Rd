\name{reg_eval}
\alias{reg_eval}
\title{
Regional evaluation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
reg_eval(species.ps, ref_area.prj, focal_list, Area_sp_in_ref, Area_ref, thres = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{species.ps}{
%%     ~~Describe \code{species.ps} here~~
}
  \item{ref_area.prj}{
%%     ~~Describe \code{ref_area.prj} here~~
}
  \item{focal_list}{
%%     ~~Describe \code{focal_list} here~~
}
  \item{Area_sp_in_ref}{
%%     ~~Describe \code{Area_sp_in_ref} here~~
}
  \item{Area_ref}{
%%     ~~Describe \code{Area_ref} here~~
}
  \item{thres}{
%%     ~~Describe \code{thres} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (species.ps, ref_area.prj, focal_list, Area_sp_in_ref, 
    Area_ref, thres = 2) 
{
    require(PBSmapping)
    Area_focal <- matrix("NA", nrow = 1, length(focal_list))
    Area_sp_in_focal <- matrix("NA", nrow = 1, length(focal_list))
    for (i in 1:length(focal_list)) {
        print(paste("calculating", focal_list[i], sep = " "))
        Area_focal[i] <- calc_area(ref_area.prj[ref_area.prj$NAME == 
            focal_list[i], ])
        try_result <- try(SpatialPolygons2PolySet(ref_area.prj[ref_area.prj$NAME == 
            focal_list[i], ]), silent = T)
        if (class(try_result)[1] == "try-error") {
            Area_sp_in_focal[i] <- "NA"
            print("reg_eval(), converting polygon to polyset error")
        }
        else {
            focal.ps <- combinePolys(SpatialPolygons2PolySet(ref_area.prj[ref_area.prj$NAME == 
                focal_list[i], ]))
            try_joinPolys <- try(joinPolys(species.ps, focal.ps, 
                "INT", maxVert = 1e+06), silent = T)
            if (class(try_joinPolys)[1] == "try-error") {
                Area_sp_in_focal[i] <- "NA"
            }
            else if (class(try_joinPolys)[1] == "NULL") {
                Area_sp_in_focal[i] <- 0
            }
            else {
                Area_sp_in_focal[i] <- sum(calcArea(combinePolys(joinPolys(species.ps, 
                  focal.ps, "INT", maxVert = 1e+06)), rollup = 1)$area)
            }
        }
    }
    DPob <- Area_sp_in_ref/Area_ref
    DPexp <- as.numeric(Area_sp_in_focal)/as.numeric(Area_focal)
    reg_ass <- DPob/DPexp
    reg_ass[which(reg_ass <= thres)] <- 1
    reg_ass[which(reg_ass > thres & reg_ass < Inf)] <- 0
    reg_ass[which(reg_ass == Inf)] <- "-9999"
    return_Matrix <- rbind(Area_sp_in_focal, as.numeric(reg_ass))
    colnames(return_Matrix) <- c(focal_list)
    return(return_Matrix)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
