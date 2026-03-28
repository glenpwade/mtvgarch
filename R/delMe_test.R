# # Initialisation: ####
# #
# setwd("C:\\Repos\\mtvgarch\\R")
# source("all_common.r")
# source("clsTVGARCH.r")
# #
#
# # Set constants
# Reps <- 20
# Tobs <- 2000
# GARCHparScale <- c(0.05,0.05,0.90)
# iterRelTol <- 1e-5  #Convergence tolerance for Iterative estimator
# estCtrl <- list(calcSE = TRUE, verbose = TRUE)
# estCtrl <- list(calcSE = FALSE, verbose = FALSE)
#
# # Check one Series ####
#
# # Set a working directory
# setwd("C:\\Repos\\LSS_LackOfIdentification")
#
# fileName = ".\\SimSourceData\\T2000_Med_VarShift.RDS"
# simData <- readRDS(fileName)
#
# # Check a few plots:
# plot(simData[,37],type='l')
#
#
# # Create the TV Specification - to be estimated later:
# st = (1:Tobs)/Tobs
# shape = tvshape$single
# # Create the TV Specification and set starting params to match the loaded Dataset
# TVspec <- MTVGARCH::tv(st,shape)
# TVspec$delta0 = runif(1,0.4,0.6)           # 0.5
# TVspec$pars["deltaN",1] = runif(1,3,5)     # 4.0
# TVspec$pars["speedN",1] = runif(1,2,3)     # 2.3
# TVspec$pars["locN1",1] = runif(1,0.4,0.6)  # 0.5
#
# TVparScale <- c(0.5,4.0,log(10),0.5)
# TVspec$optimcontrol$parscale <- TVparScale
# TVspec$optimcontrol$ndeps <- c(1e-3,1e-5,1e-5,1e-3)
# #TV@delta0free <- FALSE
#
# # 1. Do initial Estimation of g(t) assuming h(t) = 1
#
# e <- simData[,37]
# TV <- MTVGARCH::estimateTV(e,TVspec,estCtrl)
# MTVGARCH::plot(TV)
#
#
# GARCH <- MTVGARCH::garch(garchtype$general)
# GARCH$pars["omega",1] = 0.05
# GARCH$pars["alpha",1] = 0.05
# GARCH$pars["beta",1] = 0.90
# GARCH$optimcontrol$parscale <- c(0.05,0.05,0.9)
#
# #GARCH@omegafree <- FALSE  # TODO:  Causes estimateGARCH() - optim failed unexpectedly and returned NULL.
#
# # e,garchObj,estimationControl,tvObj
#
# garchObj <- GARCH
# estimationControl <- estCtrl
# tvObj <- TV
# # --
# # Use if debugging deep into .loglik()
# pars = optimpars
# # --
# GARCH <- MTVGARCH::estimateGARCH(e,GARCH,estCtrl)
# summary(GARCH)
#
# # estimateGARCH() - optim failed unexpectedly and returned NULL. Check the optim controls & starting params
#
#
#
# this$optimcontrol
