\name{find_over_num}
\alias{find_over_num}
\title{
%%  ~~function to do ... ~~
Calculate overlay numbers of two polygons
}
\description{
Find the intersection of two polygons and calculate the overlay number of items 
within background layer
}
\usage{
find_over_num(target_layer, background_layer)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{target_layer}{
    The target layer is a projected polygon
}
  \item{background_layer}{
    The background layer is a projected polygon
}
}
\details{
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
}
\author{
Lin, Cheng-Tao
}
\note{
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
function (target_layer, background_layer) 
{
    require(sp)
    intersection.lst <- sp::over(target_layer, background_layer, 
        returnList = T)
    return(length(unique(unlist(intersection.lst))))
  }
}
\keyword{ overlay }
\keyword{ intersection }
