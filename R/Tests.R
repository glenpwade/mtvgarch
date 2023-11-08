# Initialise ----
library(MTVGARCH)
# Generate test data ----
# Generate data for CCC Tests - unit tests
# Need data with standard Garch, some TV & correlation

Tobs = 2000
NrSeries = 3

st = seq(0,1,length.out=Tobs)
TV1 = tv(st,tvshape$single)
TV1$pars["locN1",1] = 0.3
#
TV2 = tv(st,tvshape$single)
TV2$pars["locN1",1] = 0.5
#
TV3 = tv(st,tvshape$single)
TV3$pars["locN1",1] = 0.5
##
G1 = garch(garchtype$general)
##
COR1 = ccc(3)
cccData <- generateRefData(NrSeries,Tobs, garchObj = G1,corrObj = COR1)

