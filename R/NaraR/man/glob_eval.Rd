\name{glob_eval}
\alias{glob_eval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
glob_eval(biogr, thres_lower = 1, thres_upper = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{biogr}{
%%     ~~Describe \code{biogr} here~~
}
  \item{thres_lower}{
%%     ~~Describe \code{thres_lower} here~~
}
  \item{thres_upper}{
%%     ~~Describe \code{thres_upper} here~~
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
function (biogr, thres_lower = 1, thres_upper = 3) 
{
    if (biogr - round(biogr) != 0 | thres_lower - round(thres_lower) != 
        0 | thres_upper - round(thres_upper) != 0) {
        print("Please input non-negative integer")
    }
    else {
        if (biogr <= thres_lower) {
            sp_distr <- c(0, "Local")
        }
        else if (biogr > thres_lower & biogr <= thres_upper) {
            sp_distr <- c(1, "Regional")
        }
        else if (biogr > thres_upper) {
            sp_distr <- c(2, "Global")
        }
        return(sp_distr)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
