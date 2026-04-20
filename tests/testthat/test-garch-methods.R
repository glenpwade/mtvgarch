#  usethis::use_test("estimateGARCH")
#

{#  Local setup for console execution ####
help("equality-expectations")

setwd("C:/Source/Repos/mtvgarch/tests/testthat")
source("setup-data.R")
library("devtools")
}

devtools::load_all()

test_that("estimateGARCH correctly handles: (e,tvObj_single,garObj,estCtrl)", {

  #plot(e,type='l')
  estCtrl <- list(calcSE=FALSE, verbose=FALSE, maxIter=1, fixStartPar=FALSE, startParAdjust=10)

  # Setup a standard tv object, with 1 transition
  Tobs <- NROW(e)
  st <- (1:Tobs) / Tobs

  # . Test 'single shape
  obj_single <- tv(st, shape = tvshape$single)

  # . Test all pars provided:
  garObj <- garch(garchtype$general)
  # Estimate object:
  set.seed(42)
  garObj <- estimateGARCH(e,obj_single,garObj,estCtrl)

  # Known pars for this data are:
  pars <- c(0.01778757, 0.08367254, 0.91563541)
  expect_equal(as.vector(garObj$Estimated$pars), pars, tolerance = 1e-6)
  #all.equal(as.vector(garObj$Estimated$pars), pars, tolerance = 1e-6)

  expect_equal(garObj$Estimated$error, FALSE)

})


test_that("estimateGARCH correctly handles: (e,garObj)", {

  #plot(e,type='l')
  # estCtrl <- list(calcSE=FALSE, verbose=FALSE, maxIter=1, fixStartPar=FALSE, startParAdjust=10)

  # Estimate object:
  garObj <- garch(garchtype$general)
  set.seed(42)
  garObj <- estimateGARCH(e,garObj)

  # Known pars for this data are:
  pars <- c(0.01778757, 0.08367254, 0.91563541)
  expect_equal(as.vector(garObj$Estimated$pars), pars, tolerance = 1e-6)
  #all.equal(as.vector(garObj$Estimated$pars), pars, tolerance = 1e-6)

  expect_equal(garObj$Estimated$error, FALSE)

})


test_that("estimateGARCH correctly handles: (e,garObj,estCtrl)", {

  #plot(e,type='l')
  estCtrl <- list(calcSE=FALSE, verbose=FALSE, maxIter=1, fixStartPar=FALSE, startParAdjust=10)

  # Estimate object:
  garObj <- garch(garchtype$general)
  set.seed(42)
  garObj <- estimateGARCH(e,garObj,estCtrl)

  # Known pars for this data are:
  pars <- c(0.01778757, 0.08367254, 0.91563541)
  expect_equal(as.vector(garObj$Estimated$pars), pars, tolerance = 1e-6)
  #all.equal(as.vector(garObj$Estimated$pars), pars, tolerance = 1e-6)

  expect_equal(garObj$Estimated$error, FALSE)

})


test_that("estimateGARCH correctly handles: (e,tvObj_single,garObj)", {

  #plot(e,type='l')
  #estCtrl <- list(calcSE=FALSE, verbose=FALSE, maxIter=1, fixStartPar=FALSE, startParAdjust=10)

  # Setup a standard tv object, with 1 transition
  Tobs <- NROW(e)
  st <- (1:Tobs) / Tobs

  # . Test 'single shape
  obj_single <- tv(st, shape = tvshape$single)
  garObj <- garch(garchtype$general)

  # Estimate object:
  set.seed(42)
  garObj <- estimateGARCH(e,obj_single,garObj)

  # Known pars for this data are:
  pars <- c(0.01778757, 0.08367254, 0.91563541)
  expect_equal(as.vector(garObj$Estimated$pars), pars, tolerance = 1e-6)
  #all.equal(as.vector(garObj$Estimated$pars), pars, tolerance = 1e-6)

  expect_equal(garObj$Estimated$error, FALSE)

})

test_that("estimateGARCH correctly returns Standard Errors", {

  #plot(e,type='l')
  estCtrl <- list(calcSE=TRUE, verbose=FALSE, maxIter=1, fixStartPar=FALSE, startParAdjust=10)

  # Setup a standard tv object, with 1 transition
  Tobs <- NROW(e)
  st <- (1:Tobs) / Tobs

  # . Test 'single shape
  obj_single <- tv(st, shape = tvshape$single)
  #obj_single$pars

  # . Test all pars provided:
  garObj <- garch(garchtype$general)

  # Estimate object:
  set.seed(42)
  garObj$optimcontrol$ndeps <- c(1e-3,1e-5,1e-3)
  garObj <- estimateGARCH(e,obj_single,garObj,estCtrl)

  # Known pars for this data are:
  se <- c(5.155321639 , 4.534447107, 0.001900697)
  expect_equal(as.vector(garObj$Estimated$se) * 1000, se, tolerance = 1e-6)
  #all.equal(as.vector(garObj$Estimated$se) * 1000, se, tolerance = 1e-6)

  expect_equal(garObj$Estimated$error, FALSE)

})

testthat_tolerance()
