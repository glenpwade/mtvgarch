\name{garch_class-class}
\Rdversion{1.1}
\docType{class}
\alias{garch_class-class}
\alias{estimateGARCH,numeric,garch_class,missing-method}
\alias{getTargetValue,numeric,missing,garch_class-method}
\alias{plot,garch_class,missing-method}
\alias{summary,garch_class-method}

\title{Class \code{"garch_class"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("garch_class", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} ~~ }
    \item{\code{h}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{nr.pars}:}{Object of class \code{"integer"} ~~ }
    \item{\code{order}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{names}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{namedList}"}, directly.
Class \code{"\linkS4class{list}"}, by class "namedList", distance 2.
Class \code{"\linkS4class{vector}"}, by class "namedList", distance 3.
}
\section{Methods}{
  \describe{
    \item{estimateGARCH}{\code{signature(e = "numeric", garchObj = "garch_class", estimationControl = "missing")}: ... }
    \item{getTargetValue}{\code{signature(e = "numeric", tvObj = "missing", garchObj = "garch_class")}: ... }
    \item{plot}{\code{signature(x = "garch_class", y = "missing")}: ... }
    \item{summary}{\code{signature(object = "garch_class")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("garch_class")
}
\keyword{classes}
