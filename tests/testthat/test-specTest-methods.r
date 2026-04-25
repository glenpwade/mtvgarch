#

{#  Local setup for console execution ####
  help("equality-expectations")

  setwd("C:/Source/Repos/mtvgarch/tests/testthat")
  source("setup-data.R")
}

devtools::load_all()

## Test the constructor(s)  ####

test_that("tv specification test: test.LM.TR2(e,tvObj,testOrder=1)", {

  # Setup a standard tv object, with 1 transition
  Tobs <- NROW(e)
  st <- (1:Tobs) / Tobs

  # . Assume no transitions
  obj_single <- tv(st, shape = tvshape$delta0only)
  garObj <- garch(garchtype$general,order = 1)

  # Estimate the no-transition model
  #estCtrl <- list(calcSE=FALSE, verbose=FALSE, maxIter=1, fixStartPar=FALSE, startParAdjust=10)
  obj_single <- estimateTV(e,obj_single,garObj)
  #
  refTests <- list()
  refTests$TR2 <- test.LM.TR2(e,obj_single,testOrder = 1)
  refTests$Robust <-  test.LM.Robust(e,obj_single,testOrder = 1)
  #
  # Generate the test statistic Distribution, for comparison: (Need Ref.Data, and Ref.Tests)
  # Generate data with Garch - test will add the 'g' from our TV object
  #

  refData <- generateRefData(nr.series = 1100,nr.obs = 2000,tvObj = obj_single,garchObj = garObj)

  simControl <- list()
  simControl$saveAs <- paste("TestStatDist-",strftime(Sys.time(),format="%Y%m%d-%H%M%S",usetz = FALSE))
  simControl$numLoops <- 1100
  simControl$numCores <- parallel::detectCores() - 1
  simControl$maxTestorder <- 2

  tsDist <- testStatDist(refData,obj_single,refTests,simControl)

  # Compare: Is our testStat In or Out of the Distribution:  In => There is likely a transition...  Try estimate a new model




  expect_equal(class(tvgarObj)[1], "tvgarch_class")
  #all.equal(class(tvgarObj)[1], "tvgarch_class")

})
