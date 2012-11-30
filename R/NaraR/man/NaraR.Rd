\name{NaraR}
\alias{NaraR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
NaraR(species.prj, ref_area.prj, focal_list, Area_sp_in_ref, Area_ref, thres = 2, gthl = 1, gthu = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{species.prj}{
%%     ~~Describe \code{species.prj} here~~
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
  \item{gthl}{
%%     ~~Describe \code{gthl} here~~
}
  \item{gthu}{
%%     ~~Describe \code{gthu} here~~
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
