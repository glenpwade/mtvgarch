\name{tv_class-class}
\Rdversion{1.1}
\docType{class}
\alias{tv_class-class}
\alias{.testStatDist,tv_class,matrix,list,missing-method}
\alias{estimateTV,numeric,tv_class,missing-method}
\alias{initialize,tv_class-method}
\alias{plot,tv_class,missing-method}
\alias{summary,tv_class-method}

\title{Class \code{"tv_class"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("tv_class", ...)}.
or simply \code{var <- tv(arg1,arg2)  }
%%  ~~ describe objects here ~~
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} ~~ }
    \item{\code{st}:}{Object of class \code{"numeric"} ~~ smooth-transition vector}
    \item{\code{g}:}{Object of class \code{"numeric"} ~~ conditional variance vector}
    \item{\code{delta0free}:}{Object of class \code{"logical"} ~~ used when estimating a multiplicative mtvgarch object.  See details for more information.}
    \item{\code{nr.pars}:}{Object of class \code{"integer"} ~~ the number of parameters in the model}
    \item{\code{nr.transitions}:}{Object of class \code{"integer"} ~~ the number of transitions in the model}
    \item{\code{Tobs}:}{Object of class \code{"integer"} ~~ the total number of observations in the sample}
    \item{\code{taylor.order}:}{Object of class \code{"integer"} ~~ an integer representing the order of Taylor expansion to use.  See details for more information.}
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
    \item{.testStatDist}{\code{signature(tvObj = "tv_class", refdata = "matrix", reftests = "list", simcontrol = "missing")}: ... an internal method that generates a distribution for the LM Tests}
    \item{estimateTV}{\code{signature(e = "numeric", tvObj = "tv_class", estimationControl = "missing")}: ... estimates the provided tv object using the data 'e' and returns a new tv object containing the initial parameters and the estimated parameters }
    \item{initialize}{\code{signature(.Object = "tv_class")}: ... intenal method to initilise a tv object}
    \item{plot}{\code{signature(x = "tv_class", y = "missing")}: ... simple plot of the conditional variance}
    \item{summary}{\code{signature(object = "tv_class")}: ... summary of a tv object}
	 }
}
\references{
Anna, Timo, some papers...
%%  ~~put references to the literature/web site here~~
}
\author{
Glen Wade
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
showClass("tv_class")
}
\keyword{classes}
