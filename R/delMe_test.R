# # Initialisation: ####
#
# # Set constants
# Reps <- 20
# Tobs <- 2000
# GARCHparScale <- c(0.05,0.05,0.90)
# iterRelTol <- 1e-5  #Convergence tolerance for Iterative estimator
# estCtrl <- list(calcSE = TRUE, verbose = TRUE)
# #estCtrl <- list(calcSE = FALSE, verbose = FALSE)
#
# # Check one Series ####
#
# # Set a working directory
# fileName = ".\\SimSourceData\\T2000_Med_VarShift.RDS"
# setwd("C:\\Repos\\LSS_LackOfIdentification")
#
# setwd("C:\\Source\\Repos\\LSS_LackOfIdentification")
#
# simData <- readRDS(fileName)
#
# # Check a few plots:
# plot(simData[,20],type='l')
#
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
# #TV@delta0free <- FALSE
#
# # 1. Do initial Estimation of g(t) assuming h(t) = 1
#
# e <- simData[,20]
# TV <- estimateTV(e,TVspec,estCtrl)
# summary(TV)
# plot(TV)
#
#
# GARCH <- garch(garchtype$general)
# GARCH$pars["omega",1] = 0.05
# GARCH$pars["alpha",1] = 0.05
# GARCH$pars["beta",1] = 0.90
# GARCH$optimcontrol$parscale <- c(0.05,0.05,0.9)
#
# #GARCH@omegafree <- FALSE  # TODO:  Causes estimateGARCH() - optim failed unexpectedly and returned NULL.
#
# #GARCH <- garch(garchtype$noGarch)
# GARCH <- estimateGARCH(e,GARCH,estCtrl)
# summary(GARCH)
#
# # estimateGARCH() - optim failed unexpectedly and returned NULL. Check the optim controls & starting params
#
#
# # TVGARCH Estimation:  ####
#
# # 2. Specify a multiplicitive TV GARCH model specification using the TV & GARCH specification above
# mod <- tvgarch(TV,garchType = garchtype$general)
# # 2.1 We need to set the Garch starting pars, before estimating the model:
# mod$garchpars["omega",1] = 0.05
# mod$garchpars["alpha",1] = 0.05
# mod$garchpars["beta",1] = 0.90
# mod$garchOptimcontrol$parscale <- c(0.05,0.05,0.90)
# mod$iterationReltol = 1e-08
# mod$varTarget <- TRUE
#
# # 3. Since we are only doing one - let's see what's going on & calc parameter se's.
# estCtrl <- list(calcSE = TRUE, verbose = TRUE)
#
# # 4. Run the 2-Step estimation
# mod_2s <- estimateTVGARCH_2Step(e,mod,estCtrl)
#
# # 4. Extract the estimated parameters:
# tvObj <- mod_2s$Estimated$tv
# garchpObj <- mod_2s$Estimated$garch
#
# # Save the results for reporting:  (Col:1 'EstimationMethod': Two-Step=1, Iterative=2), (Col10: 'EstimationError': 0(FALSE) / 1(TRUE))
# pars <- matrix(NA,1,10)
# colnames(pars) <- c("EstMethod","NrIterations","d0","d1","spd","loc","omega","alpha","beta","EstError")
#
# pars[1,] <- c(1,mod_2s@iterations,tvObj$delta0,tvObj$pars[1:3,1],garchpObj$pars,as.numeric(!(mod_2s$Estimated$converged)) )
# pars
#
#
#
# # 5. Run the iterative estimation
# # 5.1 But We don't need to calculate statistics for each iteration
# estCtrl <- list(calcSE = FALSE, verbose = TRUE)
#
# mod$iterationReltol = 1e-10
# mod_iter <- estimateTVGARCH_Iterate(e,mod,estCtrl)
#
# # 4. Extract the estimated parameters:
# tvObj <- mod_iter$Estimated$tv
# garchpObj <- mod_iter$Estimated$garch
# pars[1,] <- c(2,mod_iter@iterations,tvObj$delta0,tvObj$pars[1:3,1],garchpObj$pars,as.numeric(!(mod_iter$Estimated$converged)) )
#
# # Note: Any failed estimations will be identified in the last column, so the unestimated parameters can be excluded
#
# # Return:
# pars
#
#
# ##
#
# fileName = ".\\SimResults\\result_T2000_LSS_paper.RDS"
#
