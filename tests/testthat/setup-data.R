


# tests/testthat/setup-data.R
# This defines a variable available to all test-*.R files

setwd("C:/Source/Repos/mtvgarch/tests/testthat")

library(datasets)
rtns <- EuStockMarkets  # Dataset packaged with R

# Setup the data & estCtrl:
ftse_returns <- diff(log(rtns[,4]))    # 4: FTSE

#summary(log_returns)
#plot(log_returns,type='l', col="blue")


fileName = "T2000_Med_VarShift.RDS"
tstData <- readRDS(fileName)

e <- as.numeric(tstData[,1])  #TODO: Update package to handle (ts)
#plot(e,type='l')
