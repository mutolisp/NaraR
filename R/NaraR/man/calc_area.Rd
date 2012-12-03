\name{calc_area}
\alias{calc_area}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Caculate polygon area
}
\description{
Calculate polygon area use PBSmapping::calcArea()
}
\usage{
calc_area(SpatialPolygon.prj)
}
\arguments{
  \item{SpatialPolygon.prj}{
  SpatialPolygon.prj is a projected spatial polygon data 
}
}
\details{
}
\value{
}
\references{
}
\author{
Cheng-Tao Lin
School of Forestry and Resource Conservation, 
National Taiwan University
}
\note{
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{calcArea}},
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
data(world_borders)


## The function is currently defined as
function (SpatialPolygon.prj) 
{
    require(PBSmapping)
    SpatialPolygon.ps <- combinePolys(SpatialPolygons2PolySet(SpatialPolygon.prj))
    area <- sum(calcArea(SpatialPolygon.ps, rollup = 1)$area)
    return(area)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
