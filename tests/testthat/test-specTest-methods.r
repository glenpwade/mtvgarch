#

{#  Local setup for console execution ####
  help("equality-expectations")

  setwd("C:/Source/Repos/mtvgarch/tests/testthat")
  source("setup-data.R")
}

devtools::load_all()

## Test the constructor(s)  ####

test_that("tv specification test: test.LM.TR2(e,tvObj) - delta0_Only", {

  # Setup a standard tv object, with 1 transition
  Tobs <- NROW(e)
  st <- (1:Tobs) / Tobs

  # . Assume no transitions
  obj_tv0 <- tv(st, shape = tvshape$delta0only)
  # . Use best guess of Garch process (we used defaults in the dataGenerate process)
  garObj <- garch(garchtype$general,order = 1)

  # . Estimate the single-transition model
  obj_tv0 <- estimateTV(e,obj_tv0,garObj)

  #try{
   refData <- readRDS("specTest_refData_Tx0.RDS")
  # }
  # catch{
  #   refData <- generateRefData(nr.series = 1100,nr.obs = 2000,tvObj = obj_tv0,garchObj = garObj)
  #   #saveRDS(refData,"specTest_refData_Tx0.RDS")
  # }
    #
  maxTestOrd <- 3
  RefTests = list()
  for (n in  1:maxTestOrd){
    RefTests[[n]] <- list()
    RefTests[[n]] <- getTestStats(e,obj_tv0,n)
  }

  simControl <- list()
  simControl$saveAs <- paste("TestStatDist-",strftime(Sys.time(),format="%Y%m%d-%H%M%S",usetz = FALSE))
  simControl$numLoops <- 1100
  simControl$numCores <- parallel::detectCores() - 1
  simControl$maxTestorder <- 3

  # Note: testStatDist() will generate a simulated test distribution and calculate p_values as follows:
  # 1. For each column of refData, estimate the null TV model, then calculate test statistics (TR2 & Robust) for each requested TestOrder
  # 2. Then compare actual RefTest statistics with the TestStat-distributions (i.e. get p-values)
  # 3. This takes approx. 3 minutes for the TV$delta0 specification running in local-cpu parallel mode on an i7 3GHz 4-core/8-processor CPU,  Win10PC

  # Generate the test statistic Distribution, for comparison: (Need Ref.Data, and Ref.Tests)
  SIMRESULT = testStatDist(refData,obj_tv0,RefTests,simControl)
  ##
  # Print the test Results:
  testResults = data.frame()
  for(n in 1:maxTestOrd){
    testResults[n,1] = n
    testResults[n,2] = SIMRESULT[[n]]$pVal_TR2
    testResults[n,3] = SIMRESULT[[n]]$pVal_ROB
  }
  knitr::kable(testResults,'pipe',digits=3,col.names = c('Test Ord','TR2','Robust'))

  # Interpret Results:
  # testResult[order] => Represents a p-value.
  # Null model = No transition, testing for existence of a transition
  # <= 0.05 implies we reject the Null at 5%, ergo... try model with a transition

  # Which testOrder has the Lowest value?
  # testOrder 1, 2 or 3 => 1 implies a single transition is most likely, 2 => a double transition, 3 => multiple level shifts

  # TR2 Tests = 0
  expect_equal(testResults[1,2], 0)
  expect_equal(testResults[2,2], 0)
  expect_equal(testResults[3,2], 0)
  # Robust Tests = 0
  expect_equal(testResults[1,3], 0)
  expect_equal(testResults[2,3], 0)
  expect_equal(testResults[3,3], 0)

})

test_that("tv specification test: test.LM.TR2(e,tvObj) - One Transition", {

  # Setup a standard tv object, with 1 transition
  Tobs <- NROW(e)
  st <- (1:Tobs) / Tobs

  # . Assume no transitions
  obj_single <- tv(st, shape = tvshape$single)
  obj_single$delta0 <- 1.0
  obj_single$pars["deltaN",1] <- 4
  obj_single$pars["speedN",1] <- log(10)

  # . Use best guess of Garch process (we used defaults in the dataGenerate process)
  garObj <- garch(garchtype$general,order = 1)

  # . Estimate the single-transition model
  obj_single <- estimateTV(e,obj_single,garObj)

  #try{
  refData <- readRDS("specTest_refData_Tx1.RDS")
  # }
  # catch{
  #   refData <- generateRefData(nr.series = 1100,nr.obs = 2000,tvObj = obj_single,garchObj = garObj)
  #   #saveRDS(refData,"specTest_refData_Tx1.RDS")
  # }
  #
  maxTestOrd <- 3
  RefTests = list()
  for (n in  1:maxTestOrd){
    RefTests[[n]] <- list()
    RefTests[[n]] <- getTestStats(e,obj_single,n)
  }

  simControl <- list()
  simControl$saveAs <- paste("TestStatDist-",strftime(Sys.time(),format="%Y%m%d-%H%M%S",usetz = FALSE))
  simControl$numLoops <- 1100
  simControl$numCores <- parallel::detectCores() - 1
  simControl$maxTestorder <- 3

  # Note: testStatDist() will generate a simulated test distribution and calculate p_values as follows:
  # 1. For each column of refData, estimate the null TV model, then calculate test statistics (TR2 & Robust) for each requested TestOrder
  # 2. Then compare actual RefTest statistics with the TestStat-distributions (i.e. get p-values)
  # 3. This takes approx. 3 minutes for the TV$delta0 specification running in local-cpu parallel mode on an i7 3GHz 4-core/8-processor CPU,  Win10PC

  # Generate the test statistic Distribution, for comparison: (Need Ref.Data, and Ref.Tests)
  SIMRESULT = testStatDist(refData,obj_single,RefTests,simControl)
  ##
  # Print the test Results:
  testResults = data.frame()
  for(n in 1:maxTestOrd){
    testResults[n,1] = n
    testResults[n,2] = SIMRESULT[[n]]$pVal_TR2
    testResults[n,3] = SIMRESULT[[n]]$pVal_ROB
  }
  knitr::kable(testResults,'pipe',digits=4,col.names = c('Test Ord','TR2','Robust'))

  # Interpret Results:
  # testResult[order] => Represents a p-value.
  # Null model = One transition, testing for existence of a second transition
  # High pValues imply we fail to reject the Null, ergo... evidence suggests model is likely identified with a single transition

  # Which testOrder has the Lowest value?
  # testOrder 1, 2 or 3 => 1 implies a single transition is most likely, 2 => a double transition, 3 => multiple level shifts

  # TR2 Tests
  expect_equal(testResults[1,2], 0.810,tolerance = 1e-2)
  expect_equal(testResults[2,2], 0.499,tolerance = 1e-2)
  expect_equal(testResults[3,2], 0.636,tolerance = 1e-2)
  # Robust Tests
  expect_equal(testResults[1,3], 0.803,tolerance = 1e-2)
  expect_equal(testResults[2,3], 0.406,tolerance = 1e-2)
  expect_equal(testResults[3,3], 0.544,tolerance = 1e-2)

})
