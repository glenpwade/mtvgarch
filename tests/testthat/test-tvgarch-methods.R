#

{#  Local setup for console execution ####
help("equality-expectations")

setwd("C:/Source/Repos/mtvgarch/tests/testthat")
source("setup-data.R")
}

devtools::load_all()

## Test the constructor(s)  ####

test_that("tvgarch constructor (not estimated): (tvObj_single,garObj)", {

  # Setup a standard tv object, with 1 transition
  Tobs <- NROW(e)
  st <- (1:Tobs) / Tobs

  # . Test 'single shape
  obj_single <- tv(st, shape = tvshape$single)

  # . Test with general garch:
  garObj <- garch(garchtype$general)

  # . Construct object:
  tvgarObj <- tvgarch(obj_single,garObj)

  expect_equal(class(tvgarObj)[1], "tvgarch_class")
  #all.equal(class(tvgarObj)[1], "tvgarch_class")

})



## Test the estimator(s)  ####

test_that("estimateTVGARCH correctly handles: (e,tvgarObj,estCtrl), maxIter=1", {

  estCtrl <- list(calcSE=FALSE, verbose=FALSE, maxIter=1, fixStartPar=FALSE, startParAdjust=10)

  # Setup a standard tv object, with 1 transition
  Tobs <- NROW(e)
  st <- (1:Tobs) / Tobs

  # . Test 'single shape
  obj_single <- tv(st, shape = tvshape$single)

  # . Test all pars provided:
  garObj <- garch(garchtype$general)

  # . Construct object:
  tvgarObj <- tvgarch(obj_single,garObj)

  # . Estimate object:
  set.seed(42)
  tvgarObj <- estimateTVGARCH(e,tvgarObj,estCtrl)

  # Known pars for this data are:
  tvpars <- c(4.5396464, 2.2763090, 0.4989851, NA)
  expect_equal(as.vector(tvgarObj$Estimated$tv$pars), tvpars, tolerance = 1e-6)
  #all.equal(as.vector(tvgarObj$Estimated$tv$pars), tvpars, tolerance = 1e-6)

  expect_equal(tvgarObj$Estimated$error, FALSE)

})

test_that("estimateTVGARCH correctly handles: (e,tvgarObj,estCtrl), maxIter=10", {

  estCtrl <- list(calcSE=FALSE, verbose=FALSE, maxIter=10, fixStartPar=FALSE, startParAdjust=10)

  # Setup a standard tv object, with 1 transition
  Tobs <- NROW(e)
  st <- (1:Tobs) / Tobs

  # . Test 'single shape
  obj_single <- tv(st, shape = tvshape$single)

  # . Test all pars provided:
  garObj <- garch(garchtype$general)

  # . Construct object:
  tvgarObj <- tvgarch(obj_single,garObj)

  # . Estimate object:
  set.seed(42)
  tvgarObj <- estimateTVGARCH(e,tvgarObj,estCtrl)

  # Known pars for this data are:
  tvpars <- c(4.1451776, 2.3582512, 0.4758208, NA)
  expect_equal(as.vector(tvgarObj$Estimated$tv$pars), tvpars, tolerance = 1e-6)
  #all.equal(as.vector(tvgarObj$Estimated$tv$pars), tvpars, tolerance = 1e-6)

  expect_equal(tvgarObj$Estimated$error, FALSE)

  expect_equal(tvgarObj@iterations, 10)

})

test_that("estimateTVGARCH correctly handles: (e,tvgarObj,estCtrl), fixStartPar=TRUE", {

  estCtrl <- list(calcSE=FALSE, verbose=FALSE, maxIter=10, fixStartPar=TRUE, startParAdjust=10)

  # Setup a standard tv object, with 1 transition
  Tobs <- NROW(e)
  st <- (1:Tobs) / Tobs

  # . Test 'single shape
  obj_single <- tv(st, shape = tvshape$single)

  # . Test all pars provided:
  garObj <- garch(garchtype$general)

  # . Construct object:
  tvgarObj <- tvgarch(obj_single,garObj)

  # . Estimate object:
  set.seed(42)
  tvgarObj <- estimateTVGARCH(e,tvgarObj,estCtrl)

  # Known pars for this data are:
  tvpars <- c(4.1451776, 2.3582512, 0.4758208, NA)
  expect_equal(as.vector(tvgarObj$Estimated$tv$pars), tvpars, tolerance = 1e-6)
  #all.equal(as.vector(tvgarObj$Estimated$tv$pars), tvpars, tolerance = 1e-6)

  expect_equal(tvgarObj$Estimated$error, FALSE)

  expect_equal(tvgarObj@iterations, 10)

})

test_that("estimateTVGARCH correctly handles: (e,tvgarObj)", {

  if(exists("estCtrl", inherits = FALSE)) { rm(estCtrl) }
  # Note: default maxIter = 10

  # Setup a standard tv object, with 1 transition
  Tobs <- NROW(e)
  st <- (1:Tobs) / Tobs

  # . Test 'single shape
  obj_single <- tv(st, shape = tvshape$single)

  # . Test all pars provided:
  garObj <- garch(garchtype$general)
  garObj$optimcontrol$ndeps <- c(1e-3,1e-5,1e-3)

  # . Construct object:
  tvgarObj <- tvgarch(obj_single,garObj)

  # . Estimate object:
  set.seed(42)
  tvgarObj <- estimateTVGARCH(e,tvgarObj)

  # Known pars for this data are:
  tvpars <- c(4.1451776, 2.3582512, 0.4758208, NA)
  expect_equal(as.vector(tvgarObj$Estimated$tv$pars), tvpars, tolerance = 1e-6)
  #all.equal(as.vector(tvgarObj$Estimated$tv$pars), tvpars, tolerance = 1e-6)

  expect_equal(tvgarObj$Estimated$error, FALSE)

  expect_equal(tvgarObj@iterations, 10)

})

test_that("estimateGARCH correctly returns Standard Errors", {

  #plot(e,type='l')
  estCtrl <- list(calcSE=TRUE, verbose=FALSE, maxIter=1, fixStartPar=FALSE, startParAdjust=10)

  # Setup a standard tv object, with 1 transition
  Tobs <- NROW(e)
  st <- (1:Tobs) / Tobs

  # . Test 'single shape
  obj_single <- tv(st, shape = tvshape$single)

  # . Test all pars provided:
  garObj <- garch(garchtype$general)
  garObj$optimcontrol$ndeps <- c(1e-3,1e-5,1e-3)

  # . Construct object:
  tvgarObj <- tvgarch(obj_single,garObj)

  # . Estimate object:
  set.seed(42)
  tvgarObj <- estimateTVGARCH(e,tvgarObj,estCtrl)

  # Known pars for this data are:
  se <- c(0.35076092, 0.13256250, 0.02345737, NaN)
  expect_equal(as.vector(tvgarObj$Estimated$tv$se), se, tolerance = 1e-6)
  #all.equal(as.vector(tvgarObj$Estimated$tv$se), se, tolerance = 1e-6)

  expect_equal(tvgarObj$Estimated$error, FALSE)

})

