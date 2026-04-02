# # TODO: Write real unit tests with test framework
#
# # # Initialise ----
# # library(MTVGARCH)
# # # Generate test data ----
# # # Generate data for CCC Tests - unit tests
# # # Need data with standard Garch, some TV & correlation
# #
# # Tobs = 2000
# # NrSeries = 3
# #
# # st = seq(0,1,length.out=Tobs)
# # TV1 = tv(st,tvshape$single)
# # TV1$pars["locN1",1] = 0.3
# # #
# # TV2 = tv(st,tvshape$single)
# # TV2$pars["locN1",1] = 0.5
# # #
# # TV3 = tv(st,tvshape$single)
# # TV3$pars["locN1",1] = 0.5
# # ##
# # G1 = garch(garchtype$general)
# # ##
# # COR1 = ccc(3)
# # cccData <- generateRefData(NrSeries,Tobs, garchObj = G1,corrObj = COR1)
# #
#
# library(MTVGARCH)
# library(knitr)
# library(foreach)
# 
#
# # Set a working directory
# setwd("C:\\Repos\\LSS_LackOfIdentification")
#
# # Set constants
# Reps <- 30
# Tobs <- 2000
# GARCHparScale <- c(0.05,0.05,0.90)
# iterRelTol <- 1e-5  #Convergence tolerance for Iterative estimator
# estCtrl <- list(calcSE = TRUE, verbose = TRUE)
#
# # Check on 1 Series ####
# fileName = ".\\SimSourceData\\T2000_LSS_paper.RDS"
# simData <- readRDS(fileName)
#
# # Check a few plots:
# plot(simData[,37],type='l')
# plot(simData[,77],type='l')
# plot(simData[,49],type='l')
# plot(simData[,55],type='l')
# plot(simData[,61],type='l')
#
# e <- simData[,37]
#
# # Set a seed if you want to control the starting param 'shake':
# set.seed(42)
#
# # Create the TV Specification - to be estimated later:
# st = (1:Tobs)/Tobs
# shape = tvshape$single
# # Create the TV Specification and set starting params to match the loaded Dataset
# TVspec <- tv(st,shape)
# TVspec$delta0 = runif(1,0.4,0.6)           # 0.5
# TVspec$pars["deltaN",1] = runif(1,3,5)     # 4.0
# TVspec$pars["speedN",1] = runif(1,2,3)     # 2.3
# TVspec$pars["locN1",1] = runif(1,0.4,0.6)  # 0.5
#
# TVparScale <- c(0.5,4.0,log(10),0.5)
# TVspec$optimcontrol$parscale <- TVparScale
# TVspec$optimcontrol$ndeps <- c(1e-3,1e-5,1e-5,1e-3)
#
# estCtrl <- list(calcSE = FALSE, verbose = TRUE)
#
# # 1. Do initial Estimation of g(t) assuming h(t) = 1
# TV <- estimateTV(e,TVspec,estCtrl)
# summary(TV)
# plot(TV)
#
# GARCH <- garch(garchtype$general)
# GARCH <- estimateGARCH(e,GARCH,estCtrl)
#
# # 2. Specify a multiplicitive TV GARCH model specification using the TV & GARCH specification above
# mod <- MTVGARCH::tvgarch(TV,garchType = garchtype$general)
# # 2.1 We need to set the Garch starting pars, before estimating the model:
# mod$garchpars["omega",1] = 0.1
# mod$garchpars["alpha",1] = 0.1
# mod$garchpars["beta",1] = 0.8
# mod$garchOptimcontrol$parscale <- c(0.1,0.1,0.8)
#
# # 3. Since we are only doing one - let's see what's going on & calc parameter se's.
# estCtrl <- list(calcSE = TRUE, verbose = TRUE)
#
# # 4. Run the 2-Step estimation
# mod_2s <- MTVGARCH::estimateTVGARCH_2Step(e,mod,estCtrl,autoConverge = FALSE)
#
# # 5. Run the iterative estimation
# # 5.1 But We don't need to calculate statistics for each iteration
# estCtrl <- list(calcSE = FALSE, verbose = TRUE)
# mod_iter <- MTVGARCH::estimateTVGARCH_Iterate(e,mod,estCtrl,autoConverge = TRUE)
#
#
#
#
