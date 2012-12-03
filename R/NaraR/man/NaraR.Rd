\name{NaraR}
\alias{NaraR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
calculate national responsibility of biodiversity
}
\description{
Calculate national responsibility of biodiversity
}
\usage{
NaraR(species.prj, ref_area.prj, focal_list, Area_sp_in_ref, Area_ref, thres = 2, gthl = 1, gthu = 3)
}
\arguments{
  \item{species.prj}{
  species is a projected polygon for assessment
}
  \item{ref_area.prj}{
  reference area is a projected polygon 
}
  \item{focal_list}{
  Focal list is a vector for target countries or regions, the name should be the same 
  as in ref_area.prj
}
  \item{Area_sp_in_ref}{
%%     ~~Describe \code{Area_sp_in_ref} here~~
}
  \item{Area_ref}{
  Total area of reference polygon (The unit should be the same, default is square kilometers)
}
  \item{thres}{
  Threshold for determining the distribution of target species in regional scale.
}
  \item{gthl}{
  The lower bound threshold of target species in global scale. 
}
  \item{gthu}{
  The upper bound threshold of target species in global scale. 
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
function (species.prj, ref_area.prj, focal_list, Area_sp_in_ref, 
    Area_ref, thres = 2, gthl = 1, gthu = 3) 
{
    species.ps <- combinePolys(SpatialPolygons2PolySet(species.prj))
    ass.glob <- glob_eval(find_over_num(species.prj, ref_bio.prj), 
        gthl, gthu)
    ass.reg <- reg_eval(species.ps, ref_area.prj, focal_list, 
        Area_sp_in_ref, Area_ref, thres = 2)
    glob_a_reg <- as.numeric(ass.glob[1]) + as.numeric(ass.reg[2, 
        ])
    glob_a_reg[glob_a_reg < 0] <- "NA"
    glob_a_reg[glob_a_reg == 3] <- "Basic"
    glob_a_reg[glob_a_reg == 2] <- "Medium"
    glob_a_reg[glob_a_reg == 1] <- "High"
    glob_a_reg[glob_a_reg == 0] <- "Very High"
    returnNaraR <- rbind(ass.reg[1, ], glob_a_reg)
    return(returnNaraR)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
