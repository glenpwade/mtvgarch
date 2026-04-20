#
# 1. Create a TVGARCH Model Spec with the parameters used in the paper.
# Note: The variance level shift is low (delta0 = 0.5, delta1=1.5) and the Persistence is low.

devtools::load_all()


# Install required Packages - once
if(FALSE){
  install.packages("MTVGARCH")  # Install using R-Studio, using the install from file option - MTVGARCH_0.9.5.7.tar.gz provided
}
library(MTVGARCH)   # Ver. 0.9.9.11

# Set a working directory
setwd("C:\\Source\\Repos\\mtvgarch")

# Gen Data ####

# Set constants
Reps <- 5
Tobs <- 2000

if(TRUE){

  st = (1:Tobs)/Tobs
  shape = tvshape$single
  # Create the tv object with default parameters
  TV <- MTVGARCH::tv(st,shape)

  # Now override parameters with desired model specification
  TV$speedopt <- speedopt$eta
  TV$delta0 = 0.5
  TV$pars["deltaN",1] = 4.0
  TV$pars["speedN",1] = log(10)
  TV$pars["locN1",1] = 0.5

  # Create the garch object with default parameters
  GARCH = MTVGARCH::garch(garchtype$general)
  # Now override parameters with desired model specification
  GARCH$pars["omega",1] = 0.05
  GARCH$pars["alpha",1] = 0.05
  GARCH$pars["beta",1] = 0.90

  ## noiseDist is a named-list describing the error-distribution and parameters
  noiseDist <- list()
  noiseDist$name = 'Normal'
  noiseDist$mean = 0
  noiseDist$sd = 1
  set.seed(1984)
  ref_Data <- generateRefData(nr.series = Reps, nr.obs = Tobs, tvObj = TV, garchObj = GARCH)

  fileName = ".\\tests\\T2000_Med_VarShift.RDS"
  saveRDS(ref_Data,fileName)

}
plot(Ref_Data[,1],type='l')
