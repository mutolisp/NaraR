\name{find_focal}
\alias{find_focal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
find_focal(ref_area.prj, focal_list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ref_area.prj}{
%%     ~~Describe \code{ref_area.prj} here~~
}
  \item{focal_list}{
%%     ~~Describe \code{focal_list} here~~
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
function (ref_area.prj, focal_list) 
{
    require(PBSmapping)
    focal_list.ps <- vector(mode = "list", length = length(focal_list))
    for (i in 1:length(focal_list)) {
        if (focal_list[i] \%in\% ref_area.prj$NAME) {
            focal_list.ps[[i]] <- combinePolys(SpatialPolygons2PolySet(ref_area.prj[ref_area.prj$NAME == 
                focal_list[i], ]))
        }
        else {
            print(paste("Error!", focal_list[i], "does not match in the reference area", 
                sep = " "))
        }
    }
    return(focal_list.ps)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
