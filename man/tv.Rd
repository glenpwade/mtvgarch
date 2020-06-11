\name{tv}
\alias{tv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
creates a new tv object of type 'tv_class'
%%  ~~function to do ... ~~
}
\description{
By default the tv object will be...
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
tv(st = "numeric", shape = "numeric")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{st}{
  the smooth transition vector
%%     ~~Describe \code{st} here~~
}
  \item{shape}{
  the shape of the transition.  Valid values are stored in the tvshape list.  Use \code{tvshape$single}, for example, to define a "single" transition.  See Details for more information.
%%     ~~Describe \code{shape} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An instance of a tv_class object is returned.
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


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")

