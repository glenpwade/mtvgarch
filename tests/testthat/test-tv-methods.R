
{#  Local setup for console execution ####
  help("equality-expectations")

  setwd("C:/Source/Repos/mtvgarch/tests/testthat")
  source("setup-data.R")
}
devtools::load_all()

## Test the constructor(s)  ####

test_that("tv() constructor correctly handles all possible shapes", {
  # Setup a standard transition variable [9]
  Tobs <- 100
  st <- (1:Tobs) / Tobs

  # 1. Test 'delta0only' shape [3]
  obj_delta0 <- tv(st, shape = tvshape$delta0only)
  expect_s4_class(obj_delta0, "tv_class")
  expect_equal(obj_delta0$shape, tvshape$delta0only)
  expect_equal(obj_delta0@nr.pars, as.integer(1))

  # 2. Test 'single' shape [3]
  obj_single <- tv(st, shape = tvshape$single)
  expect_s4_class(obj_single, "tv_class")
  expect_equal(obj_single$shape, tvshape$single)
  expect_equal(obj_single@nr.pars, as.integer(4))

  # 3. Test 'double' shape [3]
  obj_double <- tv(st, shape = tvshape$double)
  expect_s4_class(obj_double, "tv_class")
  expect_equal(obj_double$shape, tvshape$double)
  expect_equal(obj_double@nr.pars, as.integer(5))

  # 4. Test 'double1loc' shape [3]
  obj_double1loc <- tv(st, shape = tvshape$double1loc)
  expect_s4_class(obj_double1loc, "tv_class")
  expect_equal(obj_double1loc$shape, tvshape$double1loc)
  expect_equal(obj_double1loc@nr.pars, as.integer(4))
})

test_that("tv() constructor initializes default slots correctly", {
  st <- (1:10) / 10
  obj <- tv(st, shape = tvshape$single)

  # Check default numeric properties defined in the class initializer [10]
  expect_equal(obj@nr.pars, as.integer(4))
  expect_equal(obj$delta0, 1.0)
  expect_true(obj@delta0free)
  expect_equal(nrow(obj$pars), 4) # Ensure the parameter matrix is 4 rows deep [5, 10]
})


## Test the estimator(s)  ####

test_that("estimateTV correctly handles: (e,tvObj_single,garObj,estCtrl)", {

  estCtrl <- list(calcSE=FALSE, verbose=FALSE, maxIter=1, fixStartPar=FALSE, startParAdjust=10)

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
  obj_single <- estimateTV(e,obj_single,garObj,estCtrl)

  # Known pars for this data are:
  pars <- c(4.5396464, 2.2763090, 0.4989851, NA)
  expect_equal(as.vector(obj_single$Estimated$pars), pars, tolerance = 1e-6)
  #all.equal(as.vector(obj_single$Estimated$pars), pars, tolerance = 1e-6)

})

test_that("estimateTV correctly handles: (e,tvObj_single,estCtrl)", {

  estCtrl <- list(calcSE=FALSE, verbose=FALSE, maxIter=1, fixStartPar=FALSE, startParAdjust=10)

  # Setup a standard tv object, with 1 transition
  Tobs <- NROW(e)
  st <- (1:Tobs) / Tobs

  # . Test 'single shape
  obj_single <- tv(st, shape = tvshape$single)

  # Estimate object:
  set.seed(42)
  obj_single <- estimateTV(e,obj_single,estCtrl)

  # Known pars for this data are:
  pars <- c(4.5396464, 2.2763090, 0.4989851, NA)
  expect_equal(as.vector(obj_single$Estimated$pars), pars, tolerance = 1e-6)
  #all.equal(as.vector(obj_single$Estimated$pars), pars, tolerance = 1e-6)

})

test_that("estimateTV correctly handles: (e,tvObj_single)", {

  if(exists("estCtrl", inherits = FALSE)) { rm(estCtrl) }

  # Setup a standard tv object, with 1 transition
  Tobs <- NROW(e)
  st <- (1:Tobs) / Tobs

  # . Test 'single shape
  obj_single <- tv(st, shape = tvshape$single)

  # Estimate object:
  set.seed(42)
  obj_single <- estimateTV(e,obj_single)

  # Known pars for this data are:
  pars <- c(4.5396464, 2.2763090, 0.4989851, NA)
  expect_equal(as.vector(obj_single$Estimated$pars), pars, tolerance = 1e-6)
  #all.equal(as.vector(obj_single$Estimated$pars), pars, tolerance = 1e-6)

})

test_that("estimateTV correctly handles: (e,tvObj_single,garObj)", {

  if(exists("estCtrl", inherits = FALSE)) { rm(estCtrl) }

  # Setup a standard tv object, with 1 transition
  Tobs <- NROW(e)
  st <- (1:Tobs) / Tobs

  # . Test 'single shape
  obj_single <- tv(st, shape = tvshape$single)
  garObj <- garch(garchtype$general)

  # Estimate object:
  set.seed(42)
  obj_single <- estimateTV(e,obj_single,garObj)

  # Known pars for this data are:
  pars <- c(4.5396464, 2.2763090, 0.4989851, NA)
  expect_equal(as.vector(obj_single$Estimated$pars), pars, tolerance = 1e-6)
  #all.equal(as.vector(obj_single$Estimated$pars), pars, tolerance = 1e-6)

})

